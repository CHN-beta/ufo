add_rules("mode.debug", "mode.release")
set_languages("cxx23")

add_requires("pkgconfig::tbb")
add_requires("pkgconfig::yaml-cpp")

target("ufo")
    set_kind("binary")
    add_includedirs("include")
    add_files("src/ufo.cpp")
    add_packages("tbb")
    add_packages("yaml-cpp")

