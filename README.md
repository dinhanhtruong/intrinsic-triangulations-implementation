# Intrinsically aCUTE

## Running the Code

Command line arguments in order:
1. path to `.xml` scenefile, absolute or relative to working directory
2. path to output with no extension. depending on render settings either a directory or a single `.png` file with that name will be created
3. whether or not to render all the video framess ("true" or "false"). if no algorithm is specified this is ignored and a single image is rendered
4. which algorithm to run, "flip" or "refine", renders initialized mesh if none specified **[not required]**
5. how many iterations to run algorithm (maxFlips for "flip" or maxInsertions for "refine") **[only required if algo specified]**
6. minAngle for refinement **[only required if running refine]**

## Compiling a Video

After you have a folder of image frames run `make_video.sh` with these command line args:
1. directory with the images
2. output filepath (should end in `.mp4`)
3. whether or not the image directory should be deleted at the end (0 or 1)