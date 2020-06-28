#ifndef _PITCH_SHIFTER_WRAP_
#define _PITCH_SHIFTER_WRAP_

#define MAX_FRAME_LENGTH 8192

/*
	变声处理，参考博客：
	http://blogs.zynaptiq.com/bernsee/pitch-shifting-using-the-ft/

	采样率不限，采样位数为16位。
	使用方法：
	1. 定义对象 PitchShifter pitchShifter;
	2. 初始化参数 比如想要男变女 pitchShifter.Init(4, 1024, 4, 48000) 
	3. 进行变声处理  pitchShifter.ProcessData(inData, dataSize, outputData);
*/
class PitchShifter
{
public:
	PitchShifter();
	~PitchShifter();

public:

	/*
		iSemitone 半音参数 -12~12之间。
		|-12		|0|			|12|
		最粗					最尖
		iFFFrameSize 傅里叶变换 处理的帧长，一般为1024,2048等2的n次方
		iOverSample 短时傅里叶变换参数，4是比较合适的变换参数，32质量最好。
		iSampleRate 输入的音频采样率

		可多次初始化，以最后一次初始化参数为准
	*/
	void Init(int iSemitone, int iFFTFrameSize, int iOverSample, int iSampleRate);
	
	/*
		psDataIn 输入数据
		iDataSize psDataIn数据大小
		psDataOut 输出数据

		psDataIn和psDataOut可以为同一内存地址。
	*/
	int ProcessData(short* psDataIn, int iDataSize, short* psDataOut);

private:
	float	m_fInFIFO[MAX_FRAME_LENGTH];
	float	m_fOutFIFO[MAX_FRAME_LENGTH];
	float	m_fFFTworksp[2 * MAX_FRAME_LENGTH];
	float	m_fLastPhase[MAX_FRAME_LENGTH / 2 + 1];
	float	m_fSumPhase[MAX_FRAME_LENGTH / 2 + 1];
	float	m_fOutputAccum[2 * MAX_FRAME_LENGTH];
	float	m_fAnaFreq[MAX_FRAME_LENGTH];
	float	m_fAnaMagn[MAX_FRAME_LENGTH];
	float	m_fSynFreq[MAX_FRAME_LENGTH];
	float	m_fSynMagn[MAX_FRAME_LENGTH];

	float	m_fPitchShifter;//变调参数
	int		m_iFFTFrameSize;//tfft处理帧长
	int		m_iOverSample;//短时傅里叶变换参数
	int		m_iSampleRate;//输入音频的采样率

	int		m_iRover;
	float	m_f2M_PI;
};


#endif