#! /usr/bin/env perl

my ($day, $month, $year, $version, $usea4) = @ARGV;

print "\@set UPDATED $day $month $year\n";
print "\@set UPDATED-MONTH $month $year\n";
print "\@set EDITION $version\n";
print "\@set VERSION $version\n";
if ($usea4 =~ /(ON|TRUE)/i) {
    print "\@afourpaper\n";
}
