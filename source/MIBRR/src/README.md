### Dependencies

I'm using a few libraries for now, but I'm not sure how many of them will stay, 
and whatever stays will not interfere with R-side of the project.

- [Boost](https://www.boost.org)
- [JSON for Modern C++](https://github.com/nlohmann/json)
- [spdlog](https://github.com/gabime/spdlog)

Except the Boost, all the other libraries can be added to the project just by their header, 
so we should not have any portability issues. For now, I just linked them since I have them 
in my system anyway.

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

### Misc.

#### Code Style

I'm pretty sure that I don't follow your coding style, but I'll try to be consistent 
with mine and at some point we can start changing it to yours. Fow now, I'm using 
[ClangFormat](https://clang.llvm.org/docs/ClangFormat.html)  and format my code based
on [LLVM Coding Standard](http://llvm.org/docs/CodingStandards.html). ClangFormat is quite 
nice and it can format your code based on different standards, e.g., Google, Microsoft, etc.

From what I see, your style is quite close to LLVM with some minor changes. ClangFormat even 
allows you to define your style, and most likely enforce it automatically. You can based the
style on LLVM, and tweak it until you are happy. It's fun! :D

### TODO:

- [ ] Including more variable, and variable types to the JSON file
    - [ ] Find a replacement for `Rcpp::List`. Either, `std::map` or `nlohmann::json` probably.
- [ ] Initial implementation and testing of the new runGibbs
- [ ] Link the `json runGibbs(...)` or `Rcpp::List runGibbs(...)` if necessary