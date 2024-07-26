/**
 * @file		melFB.cpp
 * @date		Mar 30, 2011
 * @author  	Gary Grutzek
 * @brief		Mel filter bank
 * @details		Mel filter bank,
 * 	   			logarithmic filterbank, mel-style but without linear spacing,
 * 	   			log spacing from lowest to highest freq
 *              added switch for bark scale, knospe 2013
 * @copyright  	Copyright (c) 2012 Gary Grutzek<br>
 * 				Cologne University of Applied Sciences<br>
 * 				<br>
 * 				This program is free software: you can redistribute it and/or modify
 *				it under the terms of the GNU General Public License as published by
 *				the Free Software Foundation, either version 3 of the License, or
 * 				(at your option) any later version.<br>
 *				This program is distributed in the hope that it will be useful,
 *				but WITHOUT ANY WARRANTY; without even the implied warranty of
 *				MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *				GNU General Public License for more details.<br>
 *				You should have received a copy of the GNU General Public License
 *				along with this program. If not, see <http://www.gnu.org/licenses/>.
 */


#include "melFB.h"
#include <cmath>
#include <iostream>

#define SAMPLE_RATE 8000.0f

//#define BARK

MelFb::MelFb(int fftSize, int numBanks, int locut, int hicut)
{
    numFilterBanks = numBanks;
    
    // calc next power of 2 
    int fftOrder = (int)ceil(logf((float)fftSize)/logf(2.0f));
    fbSize = (1 << fftOrder)/2;
    
    hiCut = hicut;
    if (hiCut > SAMPLE_RATE)
        throw "HiCut Parameter too high";
    
    loCut = locut;
    fbank = NULL;
    
    if (numFilterBanks == 0) 
        throw "Number of Filterbanks is zero";
    
    calcMelFB();
    
}

MelFb::~MelFb()
{
    if (fbank != NULL)
    {
		for (int j = 0; j < numFilterBanks; j++) 
			delete [] fbank[j] ; 
		delete [] fbank; 
	} 
    //    ippsFree(fLo);
    //    ippsFree(fHi);
}


void MelFb::calcMelFB() 
{
    
    // size of frequency bins in Hz => fs/fftsize 
    Ipp32f binSize = SAMPLE_RATE/(fbSize*2.0f-1); 
    // centerfrequncy of bands
    Ipp32f *fCenter = ippsMalloc_32f(numFilterBanks);
    // lower bounds of bands
    fLo = ippsMalloc_32f(numFilterBanks);
    // higher bounds of bands
    fHi = ippsMalloc_32f(numFilterBanks);
    // mel-equidistant delta between centerfrequencies
    Ipp32f delta = (Ipp32f)(linToMel(hiCut) - linToMel(loCut)) / (numFilterBanks + 1);
    
    // multidimensional array for filterbank
    fbank = new Ipp32f *[numFilterBanks]; // alloc some pointers ...
    for (int i = 0; i < numFilterBanks ; i++) 
        fbank[i] = new Ipp32f[fbSize]; // ... and alloc an array for each pointer   
    
    // get mel equidistant center frequencies
    Ipp32f melCenter = linToMel(loCut);
    fLo[0] = roundf(loCut/binSize);    
    
    for (int i=0; i<numFilterBanks; i++) 
    {
        melCenter += delta;
        // calc center frequency and map to next fft bin
        fCenter[i] = roundf(melToLin(melCenter)/binSize);
        
        // set lower and higher bounds of neighbour bands
        if (i > 0)
            fHi[i-1] = fCenter[i];
        
        if (i < numFilterBanks-1)
            fLo[i+1] = fCenter[i];
    }
    fHi[numFilterBanks-1] = hiCut/binSize;
    
    // calc triangular filter for each band
    Ipp32f x=0;
    for (int band=0; band<numFilterBanks; band++) 
    {
        for (int k=0; k<fbSize; k++) 
        {
            if (k > fLo[band] && k <= fCenter[band]) {
                fbank[band][k] = x;
                x += 1/(fCenter[band] - fLo[band]); 
            }
            else if (k > fCenter[band] && k <= fHi[band]) {
                fbank[band][k] = x;
                x -= 1/(fHi[band] - fCenter[band]);
            }
            else 
                fbank[band][k] = 0;
        }
        x = 0;
    }
    
    // free memory   
    ippsFree(fCenter);

}


void MelFb::applyMelFB(Ipp32f **output, const Stft& stft)
{
    // the non-sparse way:
    //   IppStatus error = ippStsNoErr;
    //    Ipp32f *filter = ippsMalloc_32f(fbSize);
    //    Ipp32f *frame = ippsMalloc_32f(stft.fftLen/2);
    //    
    //    for (int window=0; window < stft.NoOfWindows; window++) {
    //        // copy window to frame-buffer
    //        ippsCopy_32f(stft.spectrogramm[window], frame, stft.fftLen/2);
    //        
    //        for (int band=0; band<numFilterBanks; band++) {
    //            ippsCopy_32f(&fbank[band][0], filter, fbSize);
    //            // multiply spectrum with filterbank
    //            ippsMul_32f_I(frame, filter, fbSize);
    //            //  cout << "band: " << band << " window: " << window << endl;
    //            // sum 
    //            ippsSum_32f(filter, fbSize, &output[window][band], ippAlgHintFast);
    //        }
    //    }
    
    // the "sparse" way is about 3-4x faster
    int maxFilterSize = (int)(ceil(fHi[numFilterBanks-1]) - floor(fLo[numFilterBanks-1]));
    Ipp32f *filterSparse = ippsMalloc_32f(maxFilterSize);
    Ipp32f *frameSparse = ippsMalloc_32f(maxFilterSize);
    
    int filterSize = 0;
    int loF = 0; 
    
    for (int window=0; window < stft.NoOfWindows; window++) {
        
        for (int band=0; band<numFilterBanks; band++) {
            
            filterSize = (int) (ceil(fHi[band]) - floor(fLo[band]));
            loF = (int)floor(fLo[band]);
            ippsZero_32f(filterSparse, maxFilterSize);
            ippsZero_32f(frameSparse, maxFilterSize);
            
            // copy window to frame-buffer
            ippsCopy_32f(&stft.spectrogramm[window][loF], frameSparse, filterSize);
            // copy "sparse" filter from lo freq to hi freq
            ippsCopy_32f(&fbank[band][loF], filterSparse, filterSize);
            // multiply spectrum with sparse filterbank
            ippsMul_32f_I(frameSparse, filterSparse, filterSize);
            // sum 
            ippsSum_32f(filterSparse, filterSize, &output[window][band], ippAlgHintFast);
        }
    }
    
    ippsFree(filterSparse);
    ippsFree(frameSparse);
    //ippsFree(filter);
    //ippsFree(frame);   
    //    if (error) 
    //        throw "applyMelFB failed";
}

#ifdef BARK
Ipp32f MelFb::linToMel(Ipp32f linFreq)
{
    // b = 26.81f/(1960+f) - 0.53
    return (26.81f*linFreq / (1960.0f+linFreq) -0.53f);
}

Ipp32f MelFb::melToLin(Ipp32f linFreq)
{
    // f = - 490(53+100b)/(-657+25*b)
    return (- 490.0f*(53.0f+100.0f*linFreq)/(-657.0f+25.0*linFreq));
}

#else
Ipp32f MelFb::linToMel(Ipp32f linFreq)
{
    // m = 2595 * log10(1+f/700)
    return (2595.0f * (logf(1.0f + linFreq / 700.0f) / logf(10.0f)));
}

Ipp32f MelFb::melToLin(Ipp32f linFreq)
{
    // f = (m/2595)*10
    return (700.0f * (powf(10.0f, (linFreq / 2595.0f)) - 1.0f));
}
#endif



