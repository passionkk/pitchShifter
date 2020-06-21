#include "Wave.h"
#include <string.h>

//////////////////////////////////////////////////////////////////////////
WaveWriter::WaveWriter()
{
	m_hWavFile = NULL;
	m_nTotalSampleSize = 0;
}

WaveWriter::~WaveWriter()
{
	Close();
}

int WaveWriter::Open(char *pFileName, int nChannel, int nSampleRate, int nBitsPerSample)
{
	if (m_hWavFile != NULL)
	{
		return -1;
	}

	//
	m_hWavFile = fopen(pFileName, "wb");
	if (m_hWavFile == NULL)
	{
		return -1;
	}

	_WAVEDESCR	descriptor;
	descriptor.riff[0] = 'R';
	descriptor.riff[1] = 'I';
	descriptor.riff[2] = 'F';
	descriptor.riff[3] = 'F';
	descriptor.size = 0;
	descriptor.wave[0] = 'W';
	descriptor.wave[1] = 'A';
	descriptor.wave[2] = 'V';
	descriptor.wave[3] = 'E';
	fwrite(&descriptor, sizeof(_WAVEDESCR), 1, m_hWavFile);

	_WAVEFORMAT fmt;
	memset(&fmt, 0, sizeof(_WAVEFORMAT));
	fmt.id[0] = 'f';
	fmt.id[1] = 'm';
	fmt.id[2] = 't';
	fmt.id[3] = ' ';
	fmt.size = 16;
	fmt.format = 0x0001;
	fmt.channels = nChannel;
	fmt.sampleRate = nSampleRate;
	fmt.byteRate   = nSampleRate*nChannel*(nBitsPerSample/8);
	fmt.blockAlign	= nChannel*nBitsPerSample/8;
	fmt.bitsPerSample = nBitsPerSample;
	fwrite(&fmt, sizeof(_WAVEFORMAT), 1, m_hWavFile);

	char id[4] = {'d', 'a', 't', 'a'};
	fwrite(id, sizeof(char), 4, m_hWavFile);
	fwrite(&descriptor.size, sizeof(int), 1, m_hWavFile);

	m_nTotalSampleSize = 0;
	
	return 0;
}

void WaveWriter::Write(char *pSampleBuf, int nSampleSize)
{
	if (m_hWavFile == NULL)
	{
		return;
	}

	fwrite(pSampleBuf, 1, nSampleSize, m_hWavFile);
	fflush(m_hWavFile);
	m_nTotalSampleSize += nSampleSize;
}

void WaveWriter::Close()
{
	if (m_hWavFile == NULL)
	{
		return;
	}

    int nTotalSize = m_nTotalSampleSize + sizeof(_WAVEDESCR)+sizeof(_WAVEFORMAT) - 8;
	fseek(m_hWavFile, 4, SEEK_SET);
    fwrite(&nTotalSize, sizeof(int), 1, m_hWavFile);
	fseek(m_hWavFile, sizeof(_WAVEDESCR) + sizeof(_WAVEFORMAT) + 4, SEEK_SET);
	fwrite(&m_nTotalSampleSize, sizeof(int), 1, m_hWavFile);
	fclose(m_hWavFile);
	m_hWavFile = NULL;

}

bool WaveWriter::IsOpen()
{
	return (m_hWavFile != NULL);
}

//////////////////////////////////////////////////////////////////////////

WaveReader::WaveReader()
: m_hWavFile(NULL)
, m_nTotalSampleSize(0)
{
    //
}


WaveReader::~WaveReader()
{
    //
}

int WaveReader::Open(char *pFileName, int &nChannel, int &nSampleRate, int &nBitsPerSample)
{
    if (m_hWavFile != NULL)
    {
        return -2;
    }

    int nRet = -1;
    do 
    {
        m_hWavFile = fopen(pFileName, "rb");
        if (m_hWavFile == NULL)
        {
            break;
        }

        _WAVEDESCR	descriptor;
        if (fread(&descriptor, 1, sizeof(_WAVEDESCR), m_hWavFile) != sizeof(_WAVEDESCR))
        {
            break;
        }

        if (
            (descriptor.riff[0] != 'R')
            || (descriptor.riff[1] != 'I')
            || (descriptor.riff[2] != 'F')
            || (descriptor.riff[3] != 'F')
            || (descriptor.wave[0] != 'W')
            || (descriptor.wave[1] != 'A')
            || (descriptor.wave[2] != 'V')
            || (descriptor.wave[3] != 'E')
            )
        {
            break;
        }

        _WAVEFORMAT fmt;
        if (fread(&fmt, 1, sizeof(_WAVEFORMAT), m_hWavFile) != sizeof(_WAVEFORMAT))
        {
            break;
        }

        if (
            (fmt.id[0] != 'f')
            || (fmt.id[1] != 'm')
            || (fmt.id[2] != 't')
            )
        {
            break;
        }

        nChannel = fmt.channels;
        nSampleRate = fmt.sampleRate;
        nBitsPerSample = fmt.bitsPerSample;

        //get until data chunk
        bool bReadDataChunkSuccess = false;
        char szChunkName[4];
        int nChunkSize;
        do 
        {
            if (fread(szChunkName, 1, 4, m_hWavFile) != 4)
            {
                break;
            }
            if (fread(&nChunkSize, 1, 4, m_hWavFile) != 4)
            {
                break;
            }
            if (szChunkName[0] == 'd' && szChunkName[1] == 'a' && szChunkName[2] == 't' && szChunkName[3] == 'a')
            {
                bReadDataChunkSuccess = true;
                m_nTotalSampleSize = nChunkSize;
                break;
            }
            else
            {
                if (nChunkSize > 24*1024*1024 || nChunkSize <= 0)
                {//24M is my idear not standard
                    break;
                }
                char *pTmp = new char[nChunkSize];
                fread(pTmp, 1, nChunkSize, m_hWavFile);
                delete []pTmp;
            }
        } while (true);

        if (!bReadDataChunkSuccess)
        {
            break;
        }

        nRet = 0;
    } while (false);

    if (nRet != 0)
    {
        if (m_hWavFile != NULL)
        {
            fclose(m_hWavFile);
            m_hWavFile = NULL;
        }
    }

    return nRet;
}

int WaveReader::Read(char *pSampleBuf, int nSampleSize)
{
    if (m_hWavFile == NULL)
    {
        return -1;
    }

    if (feof(m_hWavFile))
    {
        return -1;
    }
 
    int nRead;
    nRead = fread(pSampleBuf, 1, nSampleSize, m_hWavFile);

    return nRead;
}

void WaveReader::Close()
{
    if (m_hWavFile != NULL)
    {
        fclose(m_hWavFile);
        m_hWavFile = NULL;
    }
}

bool WaveReader::IsOpen()
{
    return (m_hWavFile != NULL);
}

bool WaveReader::IsReadEnd()
{
    if (m_hWavFile != NULL)
    {
        return feof(m_hWavFile);
    }

    return false;
}

