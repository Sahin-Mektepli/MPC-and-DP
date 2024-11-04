# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/build/dependencies/cryptopp-cmake/cryptopp"
  "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/build/_deps/cryptopp-build"
  "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/build/_deps/cryptopp-subbuild/cryptopp-populate-prefix"
  "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/build/_deps/cryptopp-subbuild/cryptopp-populate-prefix/tmp"
  "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/build/_deps/cryptopp-subbuild/cryptopp-populate-prefix/src/cryptopp-populate-stamp"
  "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/build/_deps/cryptopp-subbuild/cryptopp-populate-prefix/src"
  "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/build/_deps/cryptopp-subbuild/cryptopp-populate-prefix/src/cryptopp-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/build/_deps/cryptopp-subbuild/cryptopp-populate-prefix/src/cryptopp-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/build/_deps/cryptopp-subbuild/cryptopp-populate-prefix/src/cryptopp-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
