#!/usr/bin/perl

#Compile GetCleft For MacOS10.4+

$CXX = "g++-4.0";

$DEFS = "";
$INCLUDES = "-I.";
$SYS = "-isysroot /Developer/SDKs/MacOSX10.4u.sdk";
$MACHOPT = "-DMACOSX_DEPLOYMENT_TARGET=10.4 -mmacosx-version-min=10.4";

$CXXFLAGS = "-g -Wall $SYS -arch i386 $MACHOPT";


$com = "$CXX Get_Cleft.c $INCLUDES $DEFS $CXXFLAGS -o GetCleft_MAC32";

print "$com\n";
`$com`;

