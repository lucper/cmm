# Algorithms for CMM

Prepare build:
```bash
# 32 bits
cmake -B build -DUSE_32=ON -DCMAKE_BUILD_TYPE=Debug
cmake -B build -DUSE_32=ON -DCMAKE_BUILD_TYPE=Release

# 64 bits
cmake -B build -DUSE_64=ON -DCMAKE_BUILD_TYPE=Debug
cmake -B build -DUSE_64=ON -DCMAKE_BUILD_TYPE=Release
```

Compile:
```bash
cmake --build build
```

Run:
```bash
./build/cmm -h
```
