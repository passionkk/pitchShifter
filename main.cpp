#include <iostream>
#include <string>
using namespace std;

#include "Wave.h"
#include "PitchShifter.h"


void TestManToWoman()
{
	std::string strFilePath = "新闻联播男声30s.wav";

	std::cout << "男声变女声测试: 声调由低到高1-12" << endl;
	std::cout << "测试素材: " << strFilePath.c_str() << endl;

	int iReadUnit = 1024 * 2;
	char szReadBuf[1024 * 2];
	int iRealReadSize = 0;

	PitchShifter pitchShifter;
	int iChannel, iSampleRate, iSampleBit; 
	for (int i = 1; i < 13; i++)
	{
		WaveReader reader;
		reader.Open((char*)strFilePath.c_str(), iChannel, iSampleRate, iSampleBit);
		//如果加载成功
		if (reader.IsOpen())
		{
			int semitones = i;	// 向上移动8个半音

			pitchShifter.Init(semitones, 1024, 4, iSampleRate);

			std::string strOutput = "男变女_" + std::to_string(i) + ".wav";
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

			cout << "男变女 " << i << " 完成" << endl;
		}
	}
}

void TestWomanToMan()
{
	std::string strFilePath = "新闻联播女声30s.wav";

	std::cout << "女声变男声测试: 声调由高到低1-12" << endl;
	std::cout << "测试素材: " << strFilePath.c_str() << endl;

	int iReadUnit = 1024 * 2;
	char szReadBuf[1024 * 2];
	int iRealReadSize = 0;

	PitchShifter pitchShifter;
	int iChannel, iSampleRate, iSampleBit;
	for (int i = 1; i < 13; i++)
	{
		WaveReader reader;
		reader.Open((char*)strFilePath.c_str(), iChannel, iSampleRate, iSampleBit);
		//如果加载成功
		if (reader.IsOpen())
		{
			int semitones = i * -1;	// 向上移动8个半音
			
			pitchShifter.Init(semitones, 1024, 4, iSampleRate);

			std::string strOutput = "女变男_" + std::to_string(i) + ".wav";
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

			cout << "女变男 " << i << " 完成" << endl;
		}
	}
}

int main(int argc, char* argv[])
{
	std::cout << "Audio PichShifter Processing " << std::endl;
	PitchShifter pitchShifter;
	int iChannel, iSampleRate, iSampleBit;
	std::string strFilePath = "新闻联播男声30s.wav";
	
	int iReadUnit = 1024 * 2;
	char szReadBuf[1024 * 2];
	int iRealReadSize = 0;

	WaveReader reader;
	reader.Open((char*)strFilePath.c_str(), iChannel, iSampleRate, iSampleBit);
	//如果加载成功
	if (reader.IsOpen())
	{
		int semitones = 4;	// 向上移动8个半音
		//float pitchShift = pow(2.0f, semitones / 12.0f);    //将半音转换为因子
		//printf("pitchShift:%f\n", pitchShift);

		std::cout << "男声变女声测试: 变声参数为 " <<semitones << endl;
		std::cout << "测试素材: " << strFilePath.c_str() << endl;

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

	strFilePath = "新闻联播女声30s.wav";
	reader.Open((char*)strFilePath.c_str(), iChannel, iSampleRate, iSampleBit);
	//如果加载成功
	if (reader.IsOpen())
	{
		int semitones = -4;	// 向上移动8个半音
		//float pitchShift = pow(2.0f, semitones / 12.0f);    //将半音转换为因子
		//printf("pitchShift:%f\n", pitchShift);

		std::cout << "女声变男声测试: 变声参数为 " << semitones << endl;
		std::cout << "测试素材: " << strFilePath.c_str() << endl;

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

	//整体测试
	//TestWomanToMan();
	//TestManToWoman();

	std::cout << "按回车键退出程序" << std::endl;
	getchar();
	return 0;
}