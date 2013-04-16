#!/usr/bin/perl

#Compile GetCleft For MacOS10.4+

$CXX = "g++-4.0";

$DEFS = "";
$INCLUDES = "-I.";
$SYS = "-isysroot /Developer/SDKs/MacOSX10.4u.sdk";
$MACHOPT = "-DMACOSX_DEPLOYMENT_TARGET=10.4 -mmacosx-version-min=10.4";

$CXXFLAGS = "-Wall $SYS -arch i386 $MACHOPT -O3";


$com = "$CXX Get_Cleft.c $INCLUDES $DEFS $CXXFLAGS -o GetCleft_MacOSX4+";

print "$com\n";
`$com`;

