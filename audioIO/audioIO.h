#pragma once
#include <vector>
#include <string>

struct AudioDeviceInfo {
  std::string name;      /*!< Character string device identifier. */
  bool probed;          /*!< true if the device capabilities were successfully probed. */
  int outputChannels;   /*!< Maximum output channels supported by device. */
  int inputChannels;    /*!< Maximum input channels supported by device. */
  int duplexChannels;   /*!< Maximum simultaneous input/output channels supported by device. */
  bool isDefault;       /*!< true if this is the default output or input device. */
  std::vector<int> sampleRates; /*!< Supported sample rates (queried from list of standard rates). */

  // Default constructor.
  AudioDeviceInfo()
    :probed(false), outputChannels(0), inputChannels(0),
       duplexChannels(0), isDefault(false) {}
};

struct DeviceParams{
	int num;
	int chans;
	int fs;
	bool set;
};

class audioIO{
public:
	audioIO();
	~audioIO();

	int setcallbackparams(int, void*, void*);
	int setdevparams(DeviceParams* devpar, int devcode);
	int startdev();
	int stopdev(int close=0);

	bool params_set;
	bool started;
	std::vector<AudioDeviceInfo> devices;

	/*
	void listdevices(dev_info*);
	void changedevice(dev_info*);
	void processfile(char*);
	*/

private:
	DeviceParams* input_dev;			// Input device parameters
	DeviceParams* output_dev;			// Output device parameters
	int			buffer_size;
	int	  algosize;
	void* algofunc;	
	void* algodata;
};
