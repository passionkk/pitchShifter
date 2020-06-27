#include <iostream>
using namespace std;

#include "Wave.h"
#include "PitchShifter.h"

int main(int argc, char* argv[])
{
	std::cout << "Audio Processing " << std::endl;
	std::cout << "֧�ֽ�����ͨ��wav��ʽ�ı������." << std::endl;

	WaveReader reader;
	int iChannel, iSampleRate, iSampleBit;
	std::string strFilePath = "26_1.wav";
	reader.Open((char*)strFilePath.c_str(), iChannel, iSampleRate, iSampleBit);	
	//������سɹ�
	if (reader.IsOpen())
	{
		int iReadUnit = 8192 * 2;
		char szReadBuf[8192 * 2];
		int iRealReadSize = 0;
		float semitones = 3;	// �����ƶ�8������
		float pitchShift = pow(2.0f, semitones / 12.0f);    //������ת��Ϊ����
		printf("pitchShift:%f\n", pitchShift);

		PitchShifter pitchShifter;
		pitchShifter.Init(semitones, 1024, 4, iSampleRate);

		WaveWriter writer;
		writer.Open("pitchShift_03_1.wav", iChannel, iSampleRate, iSampleBit);
		for (; ;)
		{
			iRealReadSize = reader.Read(szReadBuf, iReadUnit);
			if (iRealReadSize > 0)
			{
				pitchShifter.DoPitchShifter((short*)szReadBuf, iRealReadSize, (short*)szReadBuf);
				writer.Write(szReadBuf, iRealReadSize);
				//Sleep(1.0*iRealReadSize / iSampleRate * 1000);
			}
			else
			{
				break;
			}
		}
		writer.Close();
	}
	std::cout << "��������˳����� \n" << std::endl;
	getchar();
	return 0;
}