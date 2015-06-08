---
layout: page
show_meta: false
title: "Source"
header:
   image_fullwidth: "RNA.jpg"
permalink: "/source/"
---

The __PROBer__ GitHub repository is [here](https://github.com/pachterlab/PROBer). Source code can also be downloaded from the [download page]({{ site.url }}/download/). Currently, __PROBer__ can be built on Linux and Mac. If building on Mac, we suggest using a package manager such as [Homebrew](http://brew.sh) to download dependencies. Homebrew is easily installed by copying and pasting the command below at a terminal prompt:

~~~
ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
~~~

Other dependencies are either included, or can be installed using package managers on the system.

#### Requirements: 

- __CMake__ version >= 3.0
    - Mac: `brew install cmake`
    - Ubuntu: `sudo apt-get install cmake`
    - CentOS: `sudo yum install cmake`
- __zlib__ (should be installed on OSX >= 10.9)
    - Mac: Should be installed by default
    - Ubuntu: `sudo apt-get install zlib1g-dev`
    - CentOS: `sudo yum install zlib-devel`

#### Download

__PROBer__ is hosted on GitHub. The source code can be obtained by cloning the repository as follows:

`git clone https://github.com/pachterlab/PROBer.git`

#### Compile

Begin by moving to the source directory:

`cd PROBer`

Make a build directory and move there:

`mkdir build`

`cd build`

Run cmake:

- `-DCMAKE_INSTALL_PREFIX:PATH=$HOME` which will put PROBer in
  `$HOME/bin` as opposed to the default (`/usr/local/bin`)

`cmake ..`

Build the code:

`make`

The __PROBer__ executable is now located in `build/src`.

To install __PROBer__  into the cmake install prefix path type:

`make install`
