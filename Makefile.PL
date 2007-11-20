use 5.008001;
use ExtUtils::MakeMaker;

use File::Basename;

WriteMakefile(
    NAME              => 'Microarray',
    VERSION_FROM      => 'lib/Microarray.pm',
    PREREQ_PM         => {  Image::ExifTool => 6.66, 
    						GD => 2.32,
    						GD::Image => 2.27,
    						Statistics::Descriptive => 2.6,
							FileHandle => 2.01,
							File::Slurp => 9999.12,
							File::Basename => 2.73,
							Test::More => 0.64,
							Test::Group => 0.07,
							Test::Differences => 0.47,
							Test::Deep => 0.096
    					},
    ($] >= 5.005 ? 
      (AUTHOR         => 'Christopher Jones <c.jones@ucl.ac.uk>') : ()),
);
