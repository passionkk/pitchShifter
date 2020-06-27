#include <iostream>
using namespace std;

#include "Wave.h"
#include "PitchShifter.h"

int main(int argc, char* argv[])
{
	std::cout << "Audio Processing " << std::endl;
	std::cout << "支持解析单通道wav格式的变调处理." << std::endl;

	WaveReader reader;
	int iChannel, iSampleRate, iSampleBit;
	std::string strFilePath = "26_1.wav";
	reader.Open((char*)strFilePath.c_str(), iChannel, iSampleRate, iSampleBit);	
	//如果加载成功
	if (reader.IsOpen())
	{
		int iReadUnit = 8192 * 2;
		char szReadBuf[8192 * 2];
		int iRealReadSize = 0;
		float semitones = 3;	// 向上移动8个半音
		float pitchShift = pow(2.0f, semitones / 12.0f);    //将半音转换为因子
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
	std::cout << "按任意键退出程序 \n" << std::endl;
	getchar();
	return 0;
}