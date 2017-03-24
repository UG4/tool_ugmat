# ugmat #

**ugmat** is a program to compare and analyse matrices in ug4's .mat format.

## INSTALLATION: ##
Install ugmat with ug4's package manager [ughub](https://github.com/UG4/ughub):

    ughub install ugmat

Or simply clone this repository to any place on your harddrive.


## BUILDING: ##
### As part of UG4 (e.g. after installing ugmat through 'ughub install ugmat'):###
Go to your ug4 build directory and execute:

	cmake -Dugmat=ON .
	make

The resulting ugmat executable will be written to ug4/bin


### From the source directory as standalone application (Unix): ###
- navigate to the ugmat source directory in a terminal
- create a 'build' directory:
	
		mkdir build
	
- call cmake from that directory:
	
		cd build
		cmake -DCMAKE_BUILD_TYPE=Release ..
		
- call 'make' and 'make install'. You may need root privileges for the second.

		make && make install
		
  or (on Ubuntu)
  	
		make && sudo make install


DEBUGGING:
To create a debug build, proceed as above and use

	'cmake -DCMAKE_BUILD_TYPE=Debug .'
	
instead of

	'cmake -DCMAKE_BUILD_TYPE=Release .'
