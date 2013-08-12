g++ -c -fpic alignment.cpp -I/usr/lib64/ruby/1.8/x86_64-linux
g++ -shared -o alignment.so alignment.o

