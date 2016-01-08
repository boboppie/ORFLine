use strict;
use warnings;

my $term = join '|', map "\Q$_\E", split ' ', pop;
my $found;

while (<>) {
    if (/^>/) {
        $found = /$term/i ? 1 : 0;
        print if $found;
        next;
    }
    print if $found;
}
