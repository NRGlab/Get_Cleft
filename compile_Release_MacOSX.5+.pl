#!/usr/bin/perl

#Compile GetCleft For MacOS10.5+

$CXX = "g++-4.2";

$DEFS = "";
$INCLUDES = "-I.";
$SYS = "-isysroot /Developer/SDKs/MacOSX10.5.sdk";
$MACHOPT = "-DMACOSX_DEPLOYMENT_TARGET=10.5 -mmacosx-version-min=10.5";

$CXXFLAGS = "-Wall $SYS -arch i386 $MACHOPT -arch x86_64 $MACHOPT -O3";


$com = "$CXX Get_Cleft.c $INCLUDES $DEFS $CXXFLAGS -o GetCleft_MacOSX5+";

print "$com\n";
`$com`;

