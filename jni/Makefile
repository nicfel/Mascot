# Copyright (C) 2015 Tomas Flouri
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact: Tomas Flouri <Tomas.Flouri@h-its.org>,
# Heidelberg Institute for Theoretical Studies,
# Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany

# Profiling
#PROFILING=-g -pg
#PROFILING=-g

# Compiler warnings
WARN=-Wall -Wsign-compare -mavx2

CPPFLAGS = -O3 -Ofast -ftree-vectorize -ffast-math -mavx2
CC = g++
CFLAGS = -Wall -O3 -Ofast -ffast-math $(WARN) -fpermissive -fPIC -I../../ -I/opt/java/include -I/opt/java/include/linux -I/Library/Java/JavaVirtualMachines/jdk1.8.0_73.jdk/Contents/Home/include -I/Library/Java/JavaVirtualMachines/jdk1.8.0_73.jdk/Contents/Home/include/darwin
LINKFLAGS = -L. $(PROFILING)
LIBS=-lm
PROG=mascot
#beast_mascot_ode_Euler2ndOrderNative.o
OBJS=	beast_mascot_distribution_MascotNative2.o \
	Mascot.o \
	Euler2ndOrder.o Euler2ndOrderCPU.o \
	Euler2ndOrderCPU2.o\
	Euler2ndOrderCPU3.o\
	Euler2ndOrderCPU4.o\
	Euler2ndOrderCPU5.o\
	Euler2ndOrderCPU6.o\
	Euler2ndOrderCPU7.o\
	Euler2ndOrderCPU8.o\
	Euler2ndOrderCPU9.o\
	Euler2ndOrderCPU10.o\
	Euler2ndOrderCPU11.o\
	Euler2ndOrderCPU12.o\
	Euler2ndOrderCPU13.o\
	Euler2ndOrderCPU14.o\
	Euler2ndOrderCPU15.o

osx: $(OBJS)
		gcc -lc++ -lm -Wl -Ofast -shared -o libmascot.jnilib $(OBJS)

windows: $(OBJS)
		gcc -lc++ -lm -Wl -Ofast -shared -o libmascot.dll $(OBJS)

linux: $(OBJS)
		g++ -lc -lstdc++ -lm -Ofast -shared -o libmascot.so $(OBJS)

$(PROG): $(OBJS)
	$(CC) $(LINKFLAGS) -o $@ $+ $(LIBS)

%.o: %.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *~ $(OBJS) gmon.out $(PROG)


