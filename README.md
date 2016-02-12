# README #

## WormPixel(formerly OneWorm) ##

### What is this repository for? ###

* WormPixel is a library of image analysis modules that can be plugged together to create custom image analysis applications 

* WormPixel was originally developed to do high throughput biological assays on the genetic model organism C. elegans ("The Worm") but can be used for any image analysis task

* The advantage that WormPixel offers over more traditional image analysis packages is that the modules are built to cluster the granular image analysis tasks into a four [useful blocks](https://sites.google.com/site/gunnarkleemann/home/wormview-package) ((1) tool calibration, (2) particle collection, (3) shape analysis, and (4) data compilation/ figure generation), thus allowing the rapid construction of image analysis application.

* Because OneWorm is open source it allows a user to access and modify the granular image analysis tasks if needed.

* OneWorm is implemented in 

* Version 1

### How do I get set up? ###

* The modules are written as Matlab functions 

* Refer to the [architecture diagram](https://sites.google.com/site/gunnarkleemann/home/wormview-package)

* Clone the repository and add it to your search path

* The output of each module is designed to take the input the previous module

* Dependencies: Matlab image processing package

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* gunnarkl@gmail.com