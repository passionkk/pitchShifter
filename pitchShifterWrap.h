#ifndef _PITCH_SHIFTER_WRAP_
#define _PITCH_SHIFTER_WRAP_


class PitchShifter
{
public:
	PitchShifter();
	~PitchShifter();

public:
	void Init();
	/*
		fPitchShift 变声主要参数
		iSampleSum 总的采样数
		iFFFrameSize 快速傅里叶 处理的帧长，一般为1024,2048等2的n次方
		oversampling factor
	*/
	int DoPitchShifter(float fPitchShift, int iSampleSum, int iFFTFrameSize, long osamp, float sampleRate, short* indata, short *outdata);
	void Uninit();
};


#endif