
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>	// includes sin
#include <string>
#include <fstream>
#include <iostream>
#include <cstring>
#include <algorithm>

// CONSTANTS ******************************

#define PI              	3.14159265358979

// Frequency of tone to be created (A = 440 Hz)
#define TONE_FREQUENCY		440

// Duration of tone to be created (3 seconds for now)
#define SECONDS				3

//  Standard sample rate in Hz
#define SAMPLE_RATE     	44100.0

//  Standard sample size in bits
#define BITS_PER_SAMPLE		16

// Standard sample size in bytes		
#define BYTES_PER_SAMPLE	(BITS_PER_SAMPLE/8)

// Rescaling factor to convert between 16-bit shorts and doubles between -1 and 1
#define MAX_SHORT_VALUE		32768

// Number of channels
#define MONOPHONIC			1
#define STEREOPHONIC		2

// Offset of the fmt chunk in the WAV header
#define FMT_OFFSET			12

using namespace std;

double* readWavFile(int *arraySize, int *channels, char *filename);
void readWavFileHeader(int *channels, int *numSamples, FILE *inputFile);
void writeWavFileHeader(int channels, int numberSamples, double outputRate, FILE *outputFile);
void writeWavFile(double *outputArray, int outputArraySize, int channels, char *filename);
size_t fwriteIntLSB(int data, FILE *stream);
int freadIntLSB(FILE *stream);
size_t fwriteShortLSB(short int data, FILE *stream);
short int freadShortLSB(FILE *stream);
double* fft(double *inputarr, int size, int dir);
void multiply(const double *array1, const double *array2, double *outputArray, int arraySize);


double* readWavFile(int *arraySize, int *channels, char *filename) {
    double *array;
  
    FILE *inputFileStream = fopen(filename, "rb");
    if (inputFileStream == NULL) {
        printf("File %s could not be opened for reading\n", filename);
        exit(-1);
    }

    int numSamples;

    readWavFileHeader(channels, &numSamples, inputFileStream); 

    printf("Channels: %d\n", *channels);
    printf("Number of samples: %d\n", numSamples);

    //check that the number of samples is greater than zero (to avoid problems with arrays)
    if (numSamples <= 0) {
        printf("The file %s doesn't contain any samples. Exiting the program.\n", filename);
        exit(0);
    }
  
    *arraySize = numSamples * (*channels);
    
    array = new double[*arraySize];
    short *intArray = new short[*arraySize];

    int count = fread(intArray, BYTES_PER_SAMPLE, numSamples, inputFileStream);

    //  for (int i=0; i<numSamples; i++) {
    //	array[i] = ((double) intArray[i])/ MAX_SHORT_VALUE;
    //  }

    int largest = 0;
    for (int i=0; i< *arraySize; i++) {
        if (abs(intArray[i]) > largest) {
            largest = abs(intArray[i]);
        }
    }

    for (int i=0; i<*arraySize; i++) {
        array[i] = ((double) intArray[i])/ largest;
    }
        


    //clean up memory
    delete [] intArray;
    //int fclose(FILE *inputFileStream);
    
    return array;
}


void readWavFileHeader(int *channels, int *numSamples, FILE *inputFile) {
    int sampleRate;
    int bytesPerSecond;
    int dataChunkSize;
  
    unsigned char buffer[64];
    fread(buffer, sizeof(unsigned char), FMT_OFFSET, inputFile);
 
    freadIntLSB(inputFile);
    int fmtSize = freadIntLSB(inputFile);
    freadShortLSB(inputFile); //audio fmt = 1
 
    *channels = freadShortLSB(inputFile);
    sampleRate = freadIntLSB(inputFile);
    bytesPerSecond = freadIntLSB(inputFile);

    int frameSize = freadShortLSB(inputFile);
    int bitRate = freadShortLSB(inputFile);
 
    //confirm bit rate is 16
    if (bitRate != BITS_PER_SAMPLE) {
        printf("Error: bit rate of provided WAV file is not 16. Exiting.");
        exit(-1);
    }

    //confirm sample rate is 44.1 kHz
    if (sampleRate != SAMPLE_RATE) {
        printf("Error: sample rate of provided WAV file is not 44.1 Hz. Exiting.");
        exit(-1);
    }

    //skip forward to dataChunkSize
    fread(buffer, sizeof(unsigned char), fmtSize - 12, inputFile);
    
    dataChunkSize = freadIntLSB(inputFile);
    printf("Data chunk size: %d\n", dataChunkSize);

    *numSamples = dataChunkSize/ (BYTES_PER_SAMPLE * (*channels));
 
}



/*
Writes the header for a WAV file with the given attributes to 
 the provided filestream
*/

void writeWavFileHeader(int channels, int numberSamples, double outputRate, FILE *outputFile) {
    // Note: channels is not currently used. You will need to add this functionality
	// yourself for the bonus part of the assignment
	
	/*  Calculate the total number of bytes for the data chunk  */
    int dataChunkSize = numberSamples * BYTES_PER_SAMPLE;
	
    /*  Calculate the total number of bytes for the form size  */
    int formSize = 36 + dataChunkSize;
	
    /*  Calculate the total number of bytes per frame  */
    short int frameSize = channels * BYTES_PER_SAMPLE;
	
    /*  Calculate the byte rate  */
    int bytesPerSecond = (int)ceil(outputRate * frameSize);

    /*  Write header to file  */
    /*  Form container identifier  */
    fputs("RIFF", outputFile);
      
    /*  Form size  */
    fwriteIntLSB(formSize, outputFile);
      
    /*  Form container type  */
    fputs("WAVE", outputFile);

    /*  Format chunk identifier (Note: space after 't' needed)  */
    fputs("fmt ", outputFile);
      
    /*  Format chunk size (fixed at 16 bytes)  */
    fwriteIntLSB(16, outputFile);

    /*  Compression code:  1 = PCM  */
    fwriteShortLSB(1, outputFile);

    /*  Number of channels  */
    fwriteShortLSB((short) channels, outputFile);

    /*  Output Sample Rate  */
    fwriteIntLSB((int)outputRate, outputFile);

    /*  Bytes per second  */
    fwriteIntLSB(bytesPerSecond, outputFile);

    /*  Block alignment (frame size)  */
    fwriteShortLSB(frameSize, outputFile);

    /*  Bits per sample  */
    fwriteShortLSB(BITS_PER_SAMPLE, outputFile);

    /*  Sound Data chunk identifier  */
    fputs("data", outputFile);

    /*  Chunk size  */
    fwriteIntLSB(dataChunkSize, outputFile);
}


/*
Creates a WAV file with the contents of the provided outputArray as the samples, and writes
it to the given filename
 */

void writeWavFile(double *outputArray, int outputArraySize, int channels, char *filename) {
    // Note: channels is not currently used. You will need to add this functionality
	// yourself for the bonus part of the assignment

  //open a binary output file stream for writing
    FILE *outputFileStream = fopen(filename, "wb");
    if (outputFileStream == NULL) {
      printf("File %s cannot be opened for writing\n", filename);
        return;
    }

    //create an 16-bit integer array to hold rescaled samples
    short *intArray = new short[outputArraySize];

    //find the largest entry and uses that to rescale all other
    // doubles to be in the range (-1, 1) to prevent 16-bit integer overflow
    double largestDouble = 1;
    for (int i=0; i< outputArraySize; i++) {
		if (abs(outputArray[i]) > largestDouble) {
			largestDouble = abs(outputArray[i]);
		}
    }

    for (int i=0; i<outputArraySize; i++) {
		intArray[i] = (short) ((outputArray[i]/largestDouble)*MAX_SHORT_VALUE);
    }
	
    int numSamples = outputArraySize;

	// actual file writing
    writeWavFileHeader(channels, numSamples, SAMPLE_RATE, outputFileStream);
    fwrite(intArray, sizeof(short), outputArraySize, outputFileStream);
    
    //clear memory from heap
    delete[] intArray;
}


//writes an integer to the provided stream in little-endian form
size_t fwriteIntLSB(int data, FILE *stream) {
    unsigned char array[4];

    array[3] = (unsigned char)((data >> 24) & 0xFF);
    array[2] = (unsigned char)((data >> 16) & 0xFF);
    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 4, stream);
}


//reads an integer from the provided stream in little-endian form
int freadIntLSB(FILE *stream) {
    unsigned char array[4];

    fread(array, sizeof(unsigned char), 4, stream);
    
    int data;
    data = array[0] | (array[1] << 8) | (array[2] << 16) | (array[3] << 24);

    return data;
}


//writes a short integer to the provided stream in little-endian form
size_t fwriteShortLSB(short int data, FILE *stream) {
    unsigned char array[2];

    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 2, stream);
}


//reads a short integer from the provided stream in little-endian form
short int freadShortLSB(FILE *stream) {
    unsigned char array[2];

    fread(array, sizeof(unsigned char), 2, stream);
    
    int data;
    data = array[0] | (array[1] << 8);

    return data;
}

int main(int argc, char **argv) {
	char *outputFileName;
    char *impulseFileName;
    char *inputFileName;
	if (argc <= 3) {
		printf("Please enter input file, impulse file and output file\n");
		exit(-1);
	}
    inputFileName = argv[1];
    impulseFileName = argv[2];
	outputFileName = argv[3];
	

    printf("Reading input array from %s....\n", inputFileName);
    int inputChannels;
    double *inputArray;
    int inputArraySize;
    
    inputArray = readWavFile(&inputArraySize, &inputChannels, inputFileName);
    cout << "input channels: " << inputChannels << endl;

    cout << "Reading input array from " << impulseFileName << endl;
    int impulseChannels;
    double *impulseArray;
    int impulseArraySize;

    impulseArray = readWavFile(&impulseArraySize, &impulseChannels, impulseFileName);
    cout << "impulse channels: " << impulseChannels << endl;

    // mono 
    int outputChannels = 1;
    int outputArraySize;

    // P = n+m-1
    outputArraySize = (inputArraySize) + (impulseArraySize) - 1;
    double *outputArray = new double[outputArraySize];


    double power = log2(outputArraySize);
    power = ceil(power);

    cout << "power " << power << endl;

    int freqSize =  pow(2,power) * 2;

    cout << "Freq size " << freqSize << endl;

    double *inputFreq = new double[freqSize];
    double *impulseFreq = new double[freqSize];

    for(int i = 0; i < freqSize; i++){  
        inputFreq[i] = 0.0;
    }

    for(int i = 0; i < freqSize; i++){
        impulseFreq[i] = 0.0;
    }
    
    int counter;
    // OPTIMIZATION 2: Partial unrolling
    for(counter = 0; counter < inputArraySize-1; counter+=2){
        inputFreq[counter*2] = inputArray[counter];
        inputFreq[2*(counter+1)] = inputArray[counter+1];
    }

    if(counter == inputArraySize-1){
        inputFreq[counter*2] = inputArray[counter];
    }


    // OPTIMIZATION 2: Partial unrolling
    for(counter = 0; counter < impulseArraySize-1; counter+=2){
        impulseFreq[counter*2] = impulseArray[counter];
        impulseFreq[2*(counter+1)] = impulseArray[counter+1];
    }

    if(counter == impulseArraySize-1){
        impulseFreq[counter*2] = impulseArray[counter];
    }

    // clear the outputArray
    for(int n = 0; n < outputArraySize; n++){
        outputArray[n] = 0.0;
    }

    double *input_fft = fft(inputFreq,freqSize,1);
    double *impulse_fft = fft(impulseFreq,freqSize,1);

    double *outputFreqArr = new double[freqSize];

    multiply(input_fft,impulse_fft,outputFreqArr,freqSize);

    double *inverse_fft = fft(outputFreqArr,freqSize,-1);

    // OPTIMIZATION 2: Partial unrolling
    for(int i = 0; i < freqSize; i+=2){
        inverse_fft[i] = inverse_fft[i] / (freqSize/2);
        inverse_fft[i+1] = inverse_fft[i+1] / (freqSize/2);
    }

    for(int i = 0; i < outputArraySize; i++){
        outputArray[i] = inverse_fft[2*i];
    }

    writeWavFile(outputArray, outputArraySize, outputChannels, outputFileName);


    printf("Finished");

}

double* fft(double *inputarr, int size, int dir)
{
 
    double *outputarr = new double[size];

    // if input contains just one element
    if (size == 2){
        return inputarr;
    }

    // roots of unity array
    double *w = new double[size/2];
    // For storing n complex nth roots of unity
    // OPTIMIZATION 1: Minimize Work by creating a constant
    double alpha_constant = -2 * M_PI / (size/2);

    for (int i = 0; i < (size/4); i++) {
        // Changed to constant
        double alpha = i * alpha_constant;
        w[2*i] = cos(alpha);
        w[2*i+1] = dir*sin(alpha);
    }
 
    double *evenarr = new double[size/2];
    double *oddarr = new double[size/2];

    for (int i = 0; i < (size/4); i++) {
 
        // even indexed coefficients
        evenarr[2*i] = inputarr[i*4];
        evenarr[2*i+1] = inputarr[i*4+1];

        // odd indexed coefficients
        oddarr[2*i] = inputarr[i*4+2];
        oddarr[2*i+1] = inputarr[i*4+3];
    }
 
    // Recursive call for even indexed coefficients
    double *even_fft = fft(evenarr,(size/2),dir);

     // Recursive call for odd indexed coefficients
    double *odd_fft = fft(oddarr,(size/2),dir);
   
 
 
    for (int i = 0; i < (size / 4); i++) {
        double real = w[2*i] * odd_fft[2*i] - w[2*i+1] * odd_fft[2*i+1];
        double imaginary = w[2*i] * odd_fft[2*i+1] + w[2*i+1] * odd_fft[2*i];

        outputarr[(2*i)] = even_fft[2*i] + real;
        outputarr[(2*i)+1] = even_fft[(2*i)+1] + imaginary;

        outputarr[(2*i)+size/2] = even_fft[2*i] - real;
        outputarr[((2*i)+1)+size/2] = even_fft[(2*i)+1] - imaginary;

    }
    return outputarr;

}

// Code from TA
void multiply(const double *array1, const double *array2, double *outputArray, int arraySize) {
    
    // iterate over the complex entries, up to n = arraySize/2
    for (int i = 0; i < arraySize; i+=2) {
        int re = i;    // index of the real parts
        int im = i+1;    // index of the imaginary parts
        
        outputArray[re] = array1[re]*array2[re] - array1[im]*array2[im]; // compute the real part of the output entry
        outputArray[im] = array1[re]*array2[im] + array1[im]*array2[re]; // Compute the imaginary part of the output entry
    
        // 1. (a + b) * (c + d)
        // 2. a * c
        // 3. b * d
        // to get (ac - bd) -> step 2 - step 3
        // to get i(bc + ad) -> step 1 - step 2

        // a = array1[re]
        // b = array1[im]

        // c = array2[re]
        // d = array2[im]


    
    }
    
}