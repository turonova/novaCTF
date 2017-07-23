#pragma once

typedef struct
{
	int nx;	//number of columns (fastest changing in map)
	int ny;	//number of rows
	int nz;	//number of sections (slowest changing in map)

	int mode;	//0 image : signed 8-bit bytes range -128 to 127
				//1 image : 16-bit halfwords
				//2 image : 32-bit reals
				//3 transform : complex 16-bit integers
				//4 transform : complex 32-bit reals
				//6 image : unsigned 16-bit range 0 to 65535


	int nxStart;	//Starting point of sub image - column (not used in IMOD)
	int nyStart;	//Starting point of sub image - row (not used in IMOD)
	int nzStart;	//Starting point of sub image - section (not used in IMOD)

	int mx;	//Grid size X
	int my;	//Grid size Y
	int mz;	//Grid size Z

	float cellDimX;	//cell x dimension
	float cellDimY;	//cell y dimension
	float cellDimZ;	//cell z dimension

	float cellAngleX;	//cell x angle in degrees - not used in IMOD
	float cellAngleY;	//cell y angle in degrees - not used in IMOD
	float cellAngleZ;	//cell z angle in degrees - not used in IMOD

	int mapC;			//axis corresponding to columns (1,2,3 for X,Y,Z) - not used in IMOD
	int mapR;			//axis corresponding to rows (1,2,3 for X,Y,Z) - not used in IMOD
	int mapS;			//axis corresponding to sections (1,2,3 for X,Y,Z) - not used in IMOD

	float dMin;			//minimum density value
	float dMax;			//maximum density value
	float dMean;		//mean density value

	short ISPG;			//space group number 0 or 1 (default=0)

	short NSYMBT;		//number of bytes used for symmetry data (0 or 80)
	int extra;			//number of bytes in extended header

	short creatorID;	//creator ID

	char extraData[30];	//extra data - not used in IMOD

	short nint;			// Number of integers per section (Agard format) or
						// number of bytes per section (SerialEM format)
	short nreal;		//Number of reals per section (Agard format) or
                        //flags for which types of short data (SerialEM format):
                        //1 = tilt angle * 100  (2 bytes)
                        //2 = piece coordinates for montage  (6 bytes)
                        //4 = Stage position * 25    (4 bytes)
                        //8 = Magnification / 100 (2 bytes)
						//16 = Intensity * 25000  (2 bytes)
                        //32 = Exposure dose in e-/A2, a float in 4 bytes
                        //128, 512: Reserved for 4-byte items
                        //64, 256, 1024: Reserved for 2-byte items
                        //If the number of bytes implied by these flags does
                        //not add up to the value in nint, then nint and nreal
                        //are interpreted as ints and reals per section

	short   sub;
	short   zfac;

	float   min2;
	float   max2;
	float   min3;
	float   max3;
	float   min4;
	float   max4;		//extra data - not used in IMOD

	short   idtype;		//0 = mono, 1 = tilt, 2 = tilts, 3 = lina, 4 = lins
	short   lens;
	short   nd1;		//for idtype = 1, nd1 = axis (1, 2, or 3)
	short   nd2;
	short   vd1;        //vd1 = 100. * tilt increment
	short   vd2;        //vd2 = 100. * starting angle

	float tiltangles[6];	// Current angles are used to rotate a model to match a
							//new rotated image.  The three values in each set are
							//rotations about X, Y, and Z axes, applied in the order
							//Z, Y, X.
							//0,1,2 = original:  3,4,5 = current


	float originX;			//origin in X used for transforms
	float originY;			//origin in Y used for transforms
	float originZ;			//origin in Z used for transforms

	char map[4];			//character string 'MAP' to identify file type
	char machST[4];			//machine stamp
	float rms;				//rms deviation of map from mean density

	int nLabel;				//number of labels being used
	char labels[10][80];	//LABEL(20,10) 10 80-character text labels

}MRCHeader;


/*	Header extension used by IMOD and FEI
	Differs for each image and thus is to be filled in rawfile.h in getMRCProjection procedure
	Is always 128 bytes long */
typedef struct
{
	float   a_tilt;			//Alpha tilt (deg)
	float   b_tilt;			//Beta tilt (deg)
	float   x_stage;		//Stage x position (Unit=m. But if value>1, unit=???m)
	float   y_stage;		//Stage y position (Unit=m. But if value>1, unit=???m)
	float   z_stage;		//Stage z position (Unit=m. But if value>1, unit=???m)
	float   x_shift;		//Image shift x (Unit=m. But if value>1, unit=???m)
	float   y_shift;		//Image shift y (Unit=m. But if value>1, unit=???m)
	float   defocus;		//Defocus Unit=m. But if value>1, unit=???m)
	float   exp_time;		//Exposure time (s)
	float   mean_int;		//Mean value of image
	float   tilt_axis;		//Tilt axis (deg)
	float   pixel_size;		//Pixel size of image (m)
	float   magnification;  //Magnification used
	float   remainder[19];	//Not used (filling up to 128 bytes)

}ExtendedIMODHeader;

