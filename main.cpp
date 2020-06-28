#include <iostream>
#include <string>
using namespace std;

#include "Wave.h"
#include "PitchShifter.h"


void TestManToWoman()
{
	std::string strFilePath = "������������30s.wav";

	std::cout << "������Ů������: �����ɵ͵���1-12" << endl;
	std::cout << "�����ز�: " << strFilePath.c_str() << endl;

	int iReadUnit = 1024 * 2;
	char szReadBuf[1024 * 2];
	int iRealReadSize = 0;

	PitchShifter pitchShifter;
	int iChannel, iSampleRate, iSampleBit; 
	for (int i = 1; i < 13; i++)
	{
		WaveReader reader;
		reader.Open((char*)strFilePath.c_str(), iChannel, iSampleRate, iSampleBit);
		//������سɹ�
		if (reader.IsOpen())
		{
			int semitones = i;	// �����ƶ�8������

			pitchShifter.Init(semitones, 1024, 4, iSampleRate);

			std::string strOutput = "�б�Ů_" + std::to_string(i) + ".wav";
			WaveWriter writer;
			writer.Open((char*)strOutput.c_str(), iChannel, iSampleRate, iSampleBit);
			for (;;)
			{
				iRealReadSize = reader.Read(szReadBuf, iReadUnit);
				if (iRealReadSize > 0)
				{
					pitchShifter.ProcessData((short*)szReadBuf, iRealReadSize, (short*)szReadBuf);
					writer.Write(szReadBuf, iRealReadSize);
				}
				else
				{
					break;
				}
			}
			writer.Close();

			cout << "�б�Ů " << i << " ���" << endl;
		}
	}
}

void TestWomanToMan()
{
	std::string strFilePath = "��������Ů��30s.wav";

	std::cout << "Ů������������: �����ɸߵ���1-12" << endl;
	std::cout << "�����ز�: " << strFilePath.c_str() << endl;

	int iReadUnit = 1024 * 2;
	char szReadBuf[1024 * 2];
	int iRealReadSize = 0;

	PitchShifter pitchShifter;
	int iChannel, iSampleRate, iSampleBit;
	for (int i = 1; i < 13; i++)
	{
		WaveReader reader;
		reader.Open((char*)strFilePath.c_str(), iChannel, iSampleRate, iSampleBit);
		//������سɹ�
		if (reader.IsOpen())
		{
			int semitones = i * -1;	// �����ƶ�8������
			
			pitchShifter.Init(semitones, 1024, 4, iSampleRate);

			std::string strOutput = "Ů����_" + std::to_string(i) + ".wav";
			WaveWriter writer;
			writer.Open((char*)strOutput.c_str(), iChannel, iSampleRate, iSampleBit);
			for (;;)
			{
				iRealReadSize = reader.Read(szReadBuf, iReadUnit);
				if (iRealReadSize > 0)
				{
					pitchShifter.ProcessData((short*)szReadBuf, iRealReadSize, (short*)szReadBuf);
					writer.Write(szReadBuf, iRealReadSize);
				}
				else
				{
					break;
				}
			}
			writer.Close();

			cout << "Ů���� " << i << " ���" << endl;
		}
	}
}

int main(int argc, char* argv[])
{
	std::cout << "Audio PichShifter Processing " << std::endl;
	PitchShifter pitchShifter;
	int iChannel, iSampleRate, iSampleBit;
	std::string strFilePath = "������������30s.wav";
	
	int iReadUnit = 1024 * 2;
	char szReadBuf[1024 * 2];
	int iRealReadSize = 0;

	WaveReader reader;
	reader.Open((char*)strFilePath.c_str(), iChannel, iSampleRate, iSampleBit);
	//������سɹ�
	if (reader.IsOpen())
	{
		int semitones = 4;	// �����ƶ�8������
		//float pitchShift = pow(2.0f, semitones / 12.0f);    //������ת��Ϊ����
		//printf("pitchShift:%f\n", pitchShift);

		std::cout << "������Ů������: ��������Ϊ " <<semitones << endl;
		std::cout << "�����ز�: " << strFilePath.c_str() << endl;

		pitchShifter.Init(semitones, 1024, 4, iSampleRate);

		WaveWriter writer;
		writer.Open("pitchShift_4.wav", iChannel, iSampleRate, iSampleBit);
		for (; ;)
		{
			iRealReadSize = reader.Read(szReadBuf, iReadUnit);
			if (iRealReadSize > 0)
			{
				pitchShifter.ProcessData((short*)szReadBuf, iRealReadSize, (short*)szReadBuf);
				writer.Write(szReadBuf, iRealReadSize);
			}
			else
			{
				break;
			}
		}
		writer.Close();
		reader.Close();
	}

	strFilePath = "��������Ů��30s.wav";
	reader.Open((char*)strFilePath.c_str(), iChannel, iSampleRate, iSampleBit);
	//������سɹ�
	if (reader.IsOpen())
	{
		int semitones = -4;	// �����ƶ�8������
		//float pitchShift = pow(2.0f, semitones / 12.0f);    //������ת��Ϊ����
		//printf("pitchShift:%f\n", pitchShift);

		std::cout << "Ů������������: ��������Ϊ " << semitones << endl;
		std::cout << "�����ز�: " << strFilePath.c_str() << endl;

		pitchShifter.Init(semitones, 1024, 4, iSampleRate);

		WaveWriter writer;
		writer.Open("pitchShift_-4.wav", iChannel, iSampleRate, iSampleBit);
		for (;;)
		{
			iRealReadSize = reader.Read(szReadBuf, iReadUnit);
			if (iRealReadSize > 0)
			{
				pitchShifter.ProcessData((short*)szReadBuf, iRealReadSize, (short*)szReadBuf);
				writer.Write(szReadBuf, iRealReadSize);
			}
			else
			{
				break;
			}
		}
		writer.Close();
		reader.Close();
	}

	//�������
	//TestWomanToMan();
	//TestManToWoman();

	std::cout << "���س����˳�����" << std::endl;
	getchar();
	return 0;
}