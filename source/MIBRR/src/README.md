### Dependencies

I'm using a few libraries for now, but I'm not sure how many of them will stay, 
and whatever stays will not interfere with R-side of the project.

- [Boost](https://www.boost.org)
- [JSON for Modern C++](https://github.com/nlohmann/json)
- [spdlog](https://github.com/gabime/spdlog)

### Build

```shell script
cd source/MIBRR/
mkdir build && cd buid
cmake ..
make -j4
```

### Run

At some point we'll be able to do and get some sort of output from it.

```shell script
./mibrr config.json
```