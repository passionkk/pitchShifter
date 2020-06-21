#ifndef _WAVE_H_
#define _WAVE_H_

#include <stdio.h>
#include <stdlib.h>

#pragma pack(1)
typedef struct __WAVEDESCR
{
	char riff[4];
	int size;
	char wave[4];

} _WAVEDESCR, *_LPWAVEDESCR;

typedef struct __WAVEFORMAT
{
	char id[4];
	int size;
	short format;
	short channels;
	int sampleRate;
	int byteRate;
	short blockAlign;
	short bitsPerSample;

} _WAVEFORMAT, *_LPWAVEFORMAT;
#pragma pack()

class WaveWriter
{
public:
	WaveWriter();
	~WaveWriter();
private:
	FILE	*m_hWavFile;
	int		m_nTotalSampleSize;
public:
	int Open(char *pFileName, int nChannel, int nSampleRate, int nBitsPerSample);
	void Write(char *pSampleBuf, int nSampleSize);
	void Close();
	bool IsOpen();
};

class WaveReader
{
public:
    WaveReader();
    ~WaveReader();
private:
    FILE    *m_hWavFile;
    int     m_nTotalSampleSize;
public:
    int Open(char *pFileName, int &nChannel, int &nSampleRate, int &nBitsPerSample);
    int Read(char *pSampleBuf, int nSampleSize);
    void Close();
    bool IsOpen();
    bool IsReadEnd();
};

#endif //_WAVE_H_
