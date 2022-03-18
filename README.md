# su-graph-ordering
Amazingly fast graph ordering

# Compilation

```bash
cd AMD
mkdir Lib
make lib 
cd ..
cd metis-5.1.0
make config prefix=$(location of su-order)/lib/
make 
make install
cd ..
git submodule update --remote --init
cd grappolo
vim InputsOutput/loadEdgeList.cpp
```
At the top of the file, add the following line:
```c++
#define DEL_ZERO_BASED
```
Exit the file
```bash
make
cd ..
make
```
The executable is located in `bin/` with the name of the `${current_branch}_order.out$`
