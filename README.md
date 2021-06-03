# Sky Generator
This is a sky generator that calculates a physically based sky which is based on an improved version of [Nishita](https://www.scratchapixel.com/lessons/procedural-generation-virtual-worlds/simulating-sky/simulating-colors-of-the-sky) 1993 single scattering model.

### Install
* Get **Python 3.8** or a later version
* Get pip library
* Get numpy library
* Get pillow library

On Arch Linux just run the following command
```
sudo pacman -S python python-pip python-numpy python-pillow
```

### Run
To run the code, enter:
```
python main.py
```

### How it works
When finishing, the program will show the rendered sky image with the OS image visualizer.
To change the sky parameters, just change them in the `properties.py` file, there you will find the sun rotation, camera position altitude along with other settings.
