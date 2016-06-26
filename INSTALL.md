# PROBer installation

Requirements
------------
- CMake version >= 2.8.12
    - Can be installed via homebrew: `brew install cmake`
- zlib (should be installed on OSX >= 10.9)

Installation
------------

1. First clone the repository:

    ```
    git clone git@github.com:pachterlab/PROBer.git
    ```

2. Move to the source directory:

    ```
    cd PROBer
    ```

3. Make a build directory and move there:

    ```
    mkdir build
    cd build
    ```

4. Execute cmake. Set installation directory by:
    - `-DCMAKE_INSTALL_PREFIX:PATH=$HOME` which will put PROBer in
       `$HOME/bin` as opposed to the default (`/usr/local/bin`)

    ```
    cmake ..
    ```

    This is only required when one of the `CMakeLists.txt` files changes or new
    source files are introduced. It will make a new set of `Makefile`s.

5. Build the code

    ```
    make
    ```

    Optionally install into the cmake install prefix path:

    ``` make install ```

6. Run the code. You can find the executables under `build/src`. If
    you have installed PROBer, you can also find the executables under `/usr/local/bin` or `$HOME/bin`.
