#define _GNU_SOURCE
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>
#include<unistd.h>

//performs spherical transformations
int star_rotate(double starvec[4], double starcvec[4]);
//tansforms equatorial coordinates to cartesian coodinates
int star_polar2cartesian(double vec[4], double right_ascension, double declination, double r);
//tansforms cartesian coordinates to equatorial coodinates
int star_cartesian2polar(double *right_ascension, double *declination, double vec[4], int equatorial);
//used in \"millenium\" mode to fill the vector according to the right ascension axis
int star_new_north_coordinates_ra2rotate_vector(double cvec[4], double right_ascension_axis, double declination_axis);
//used in \"millenium\" mode to fill the vector according to the declination axis
int star_new_north_coordinates_dec2rotate_vector(double cvec[4], double right_ascension_axis, double declination_axis);
//performs spherical transformations for both right ascension and declination axis
int star_future_position(double cvec[4], double vec[4], double right_ascension_axis, double declination_axis);
int equation_of_time(double *EoT, int day, int month, int year, double *precession_mvt_year, double *perihelion_mvt_year, double declination_axis, int T, int calendar, double sideral_year);
//prints help
void print_help();

int main (int argc, char** argv)
{
	double starvec[4] = {0.0L};
	double starcvec[4] = {0.0L};
	double tropic_year = 365.2421898L, sideral_year = 365.256363004L, anomalistic_year = 365.259635864L;
	double cycle_precession = tropic_year / (sideral_year - tropic_year);
	double cycle_perihelion = tropic_year / (anomalistic_year - sideral_year);
	double perihelion_mvt_year = 2.0L * M_PI / cycle_perihelion, highest_obliquity = 24.5044L, lowest_obliquity = 22.0425L;
	double precession_mvt_year = 2.0L * M_PI / cycle_precession;
	double day = -1.0L, month = -1.0L, year = -1.0L;
	char *list_options = "hz:y:e:c:";
	int option, date = 0, calendar = 0;
	double right_ascension = 0.0L, right_ascension_axis = 0.0L, declination = 0.0L, declination_axis = 0.0L;
	opterr = 0;
	while((option = getopt(argc, argv, list_options)) != -1)
	{
		switch(option)
		{
			case 'h':
				print_help();
				return EXIT_SUCCESS;
			break;
			case 'y':
				sscanf(optarg, "%lf", &right_ascension);
			break;
			case 'z':
				sscanf(optarg, "%lf", &declination);
			break;
			case 'e':
				if(sscanf(optarg, "%lf %lf %lf", &day, &month, &year) != EOF)
				{
					date = 1;
				}else
				{
					printf("You have to provide all the arguments for e option\n");
					fprintf(stderr, "execute ./stars_shift -h for help\n");
					return EXIT_FAILURE;
				}
			break;
			case 'c':
				sscanf(optarg, "%d", &calendar);
			break;
			case '?':
				fprintf(stderr, "bad option: %c\n", optopt);
				fprintf(stderr, "execute ./stars_shift -h for help\n");
				return EXIT_FAILURE;
			break;
		}
	}

	if((!right_ascension && declination) || (right_ascension && !declination))
	{	
		fprintf(stderr, "argument(s) missing\n");
		fprintf(stderr, "execute ./stars_shift -h for help\n");
		return EXIT_FAILURE;
	}
	if(!date && right_ascension)
	{
		fprintf(stderr, "option -e is missing\n");
		fprintf(stderr, "execute ./stars_shift -h for help\n");

		return EXIT_FAILURE;
	}
	if(optind != argc)
	{
		fprintf(stderr, "too much arguments:\n");
		while(optind != argc)
			fprintf(stderr, "\t - %s\n", argv[optind++]);
		fprintf(stderr, "execute ./stars_shift -h for help\n");
		return EXIT_FAILURE;
	}

	// Equation of time calculus
	if(date){
		if(!calendar)
		{
			fprintf(stderr, "option -c for calendar must be specified with -e option\n");
			fprintf(stderr, "execute ./stars_shift -h for help\n");
			return EXIT_FAILURE;
		}
		
		// Gaperihelion_daya radius calculus
		right_ascension_axis = precession_mvt_year * (year - 2017.0L);
		// obliquity for the year over 41000 years... between 21.1° & 24.5°.
		double obliquity_radius = - 2.0L * M_PI / 41000.0L * (year - 2890.02L);
		declination_axis = ((highest_obliquity + lowest_obliquity) / 2.0L + (highest_obliquity - lowest_obliquity) / 2.0L  * (double)(sinl((double)(obliquity_radius)))) * M_PI / 180.0L;
	
		// Tropic year calculus
		int y = year,m = month,d = day;
		if(m < 3)
		{
			y -= 1;
			m += 12;
		}
		int JJ = (int)(365.25L*( y + 4716 ))+ (int)(30.6001L*(m + 1))+ d + 2 - (int)(y / 100) + (int)(y / 400) - 1524.5L;
		int T = (double)(JJ - 2451545) / 36525.0L;
		int t = T /10.0L;
		double tropic_year = 365.2421905166L - t * 0.000061560L - pow(t,2) * 0.0000000684L + 
			pow(t,3) * 0.0000002630L + pow(t,4) * 0.0000000032L;
		
		printf("For the day %d/%d/%d\n", (int)(day), (int)(month) , (int)(year));
		printf("\tTropic Year Value:\n");
		printf("\t\t%lf days\n", tropic_year);
		printf("\tObliquity Value:\n");
		printf("\t\t%Lf°\n", declination_axis * 180.0L / M_PI);
		
		double EoT;
		equation_of_time(&EoT, day, month, year, &precession_mvt_year, &perihelion_mvt_year, declination_axis, T, calendar, sideral_year);
		printf("\tEquation Of Time\n");
		int min = EoT * 4.0L;
		int sec = (int)(EoT * 4.0L * 60.0L) % 60;
		sec = (sec < 0)?-sec:sec;
		printf("\t\tIs by experimental method:\n");
		printf("\t\t\t%d min %d sec\n", min, sec);
		
		if(!right_ascension)
			return EXIT_SUCCESS;
	}

	right_ascension = 2.0L * M_PI - right_ascension * M_PI / 12.0L;
	declination = declination * M_PI / 180.0L;
	
	star_polar2cartesian(starvec, right_ascension, declination, 1.0L);
	
	/////////////////////////////////////////////////////////////////
	// Logodromic Debug !! Have to be coperihelion_dayented for INITIAL WORK !! 
//	printf("Equatorial logodromie before processing\n=> Must be equal to the declination after processing:\n");
//	printf("%.10lf\n",  (double)(acosl(cosl(-16.7131388888L*M_PI/180.0L)*cosl(-25.546927778L*M_PI/180.0L)*cosl(6.7525694444L*M_PI/24.0L 
//		- 2.808875L*M_PI/24.0L)+sinl(-16.7131388888L*M_PI/180.0L)*sinl(-25.546927778L*M_PI/180.0L)) * 180.0L / M_PI));
	/////////////////////////////////////////////////////////////////

	star_future_position(starcvec, starvec, right_ascension_axis, declination_axis);

	// translate the star's coordinate from cartesian to polar equatorial baseline
	star_cartesian2polar(&right_ascension, &declination, starvec, 1);

	if(right_ascension < 0)
		right_ascension = - right_ascension;
	else if(right_ascension > 0)
		right_ascension = 2.0L * M_PI - right_ascension;

	printf("Results for the star:\n");
	float hours = right_ascension * 12.0L / M_PI;
	float minutes = (hours - (float)((int)(hours))) * 60.0f;
	float seconds = (minutes - (float)((int)(minutes))) * 60.0f;
	printf("\t- Right Ascension: %dh%d'%d''\n\t- Declination: %.2Lf°\n", (int)(hours), (int)(minutes), (int)(seconds), declination * 180.0L / M_PI);
	
	return EXIT_SUCCESS;
}

void print_help()
{
	printf("\n########################################## HELP ###########################################\n");
	printf("This program converts stars's coordinates from a baseline to another\n");
	printf("in order to represent stars over the differents Ages our earth know.\n");
	printf("It also calculate the equation of time given a date ( J/M/A ),\n");
	printf("and the duration of a tropic year.\n");
	printf("The position of stars over time is the result of two cycles:\n");
	printf("\t1) the precession of equinoxes ( ~ 25770 years period ).\n");
	printf("\t2) the variations of obliquity ( ~ 41000 years period ).\n");
	printf("The third cycle is the movement of perihelion due to ellipticity ( ~ 112000 years cycle ).\n");
	printf("The obliquity, the precession & the excentricity permit to calculate the equation of time.\n");
	printf("Those cycles produce the pole shift and changes in the equation of time.\n");
	printf("BE CAREFULL: The equation of time is accurate only for the 3000 next years!!\n");
	printf("OPTIONS:\n\n");
	printf("-c : an integer that represents the calendar used:\n");
	printf("\t1: julian calendar\n");
	printf("\t2: gregorian calendar\n");
	printf("\t3: enhanced gregorian calendar that substract 1 day each 4000 & 20000 years\n");
	printf("\tto represent a tropical year of 365.2422 days.\n");
	printf("\tHave to be used when -e option is specified.\n");
	printf("-y : the right ascension of the star you want the coordinates for the Age considered.\n");
	printf("must be in fraction of hours.\n\n");
	printf("-z : the declination of the star you want the coordinates for the Age considered.\n");
	printf("must be in fraction of degrees.\n\n");
	printf("-e : a string with the day, the month and the year( separated by space ) from which you want the correction\n");
	printf("in minutes and seconds in the equation of time, or the date ( with the same format )\n");
	printf("for the new baseline from which you want your star to be represented\n");
	printf("-h : print this document\n\n");
	printf("\nExample for the 1st of January 5000:\n");
	printf("\t./stars_shift -c1 -e \"1 1 5000\" -y RIGHT_ASCENSON -z DECLINATION\n");
	printf("As long as -e option is mentioned, the equation of time is calculated\n");
	printf("We can see that Véga is our polar star by the year 14000 ( with the proper motion )\n");
	printf("\t./stars_shift -c1 -e \"1 1 13510\" -y 18.783058333 -z 39.021525\n");
	printf("\n###########################################################################################\n");
}

int equation_of_time(double *EoT, int day, int month, int year, double *precession_mvt_year, double *perihelion_mvt_year, double declination_axis, int T, int calendar, double sideral_year)
{

	double months[] = {31.0L,28.0L,31.0L,30.0L,31.0L,30.0L,31.0L,31.0L,30.0L,31.0L,30.0L,31.0L};
	double tropic_days = day;
	for(int i = 0; i < month - 1; i++)
		tropic_days += months[i];
	
	if(calendar == 0)
	{
		tropic_days -= 13.0L;
		for(int i = 2099; i < (int)(year); i++)
		{
			if(i % 100 == 0)
			{
				tropic_days -= 1.0L;
				if(i % 400 == 0)
				{
					tropic_days += 1.0L;
					
				}
			}
		}
	}
	if(calendar == 0 || calendar == 1)
	{
		for(int i = 3999; i < (int)(year); i++)
		{
			if(i % 4000 == 0)
			{
				tropic_days -= 1.0L;
				if(i % 20000 == 0)
				{
					tropic_days -= 1.0L;
				}
			}
		}
	}
	if(tropic_days < 0.0L)
	{
		year -= 1;
		tropic_days=365.0L - tropic_days;
	}else if(tropic_days > 365.0L)
	{
		year += 1;
		tropic_days=(double)((int)(tropic_days) % 365);
	}
	double perihelion_mvt = (year - 2017.0L) * (*perihelion_mvt_year);
	double precession_mvt = (*precession_mvt_year) * (year - 2017.0L);
	double mean_movement = 2.0L * M_PI / sideral_year;
	double e = 0.016708634L - 0.000042037L * (double)(T) - 0.0000001267L * (double)(pow(T,2));
	double longitude_perihelion =  (279.69668L + 36000.76892L * (double)(T) + 0.0003025L * (double)(pow(T,2))) * M_PI / 180.0L;
	double perihelion_day = mean_movement * (tropic_days - 4.0L);
	double mm = mean_movement * tropic_days;
	// Equation of Time with experimental method
	*EoT = (2.0L * e * sinl(perihelion_day - precession_mvt - perihelion_mvt) -
		pow(tanl(declination_axis / 2.0L), 2) * sinl((longitude_perihelion + mm + 
			2.0L * e * sinl(2.0L * longitude_perihelion + 3.0L * mm) - 
				2.0L * e * sinl(2.0L * longitude_perihelion + mm)) * 2.0L)) * 180.0L / M_PI;
	return EXIT_SUCCESS;

}

int star_future_position(double cvec[4], double vec[4], double right_ascension_axis, double declination_axis)
{
	// computes rotate vector from right ascension axis to represent stars according to the ecliptic
	star_new_north_coordinates_ra2rotate_vector(cvec, 0.0L, (23.0L + 26.0L / 60.0L + (12.087L + 0.4863L * 3.0L)/ 3600.0L) * M_PI / 180.0L);

	// performs transformations
	star_rotate(vec, cvec);

	// computes rotate vector from declination axis
	star_new_north_coordinates_dec2rotate_vector(cvec, right_ascension_axis, M_PI / 2.0L);
	
	// performs transformations
	star_rotate(vec, cvec);
	
	// computes rotate vector from right ascension axis
	star_new_north_coordinates_ra2rotate_vector(cvec, M_PI, declination_axis);

	// perform transformations
	star_rotate(vec, cvec);
	
	return EXIT_SUCCESS;
}

int star_new_north_coordinates_ra2rotate_vector(double cvec[4], double right_ascension_axis, double declination_axis)
{
	star_polar2cartesian(cvec, right_ascension_axis, 0.0L, 1.0L);
	cvec[0] = cosl( declination_axis / 2.0L);
	for(int i = 1; i < 4; i++)
		cvec[i] = sinl( declination_axis / 2.0L) * cvec[i];

	return EXIT_SUCCESS;
}

int star_new_north_coordinates_dec2rotate_vector(double cvec[4], double right_ascension_axis, double declination_axis)
{
	star_polar2cartesian(cvec, 0.0L, declination_axis, 1.0L);
	cvec[0] = cosl( - right_ascension_axis / 2.0L);
	for(int i = 1; i < 4; i++)
		cvec[i] = sinl( - right_ascension_axis / 2.0L) * cvec[i];

	return EXIT_SUCCESS;
}

int star_rotate(double starvec[4], double starcvec[4])
{
	double signs[4][4] = {{1.0L, -1.0L, -1.0L, -1.0L},
				{1.0L, 1.0L, -1.0L, 1.0L},
				{1.0L, 1.0L, 1.0L, -1.0L},
				{1.0L, -1.0L, 1.0L, 1.0L}};
	int values[4][4] = 	{{0, 1, 2, 3},
				{1, 0, 3, 2},
				{2, 3, 0, 1},
				{3, 2, 1, 0}};

	double startmp[4] = {0.0L};

	for(int i = 0; i < 4; i++)
		for(int j = 0; j < 4; j++)
			startmp[i] +=  starcvec[values[i][j]] * starvec[j] * signs[i][j];
	
	for(int i = 1; i < 4; i++)
		starcvec[i] = -starcvec[i];

	for(int i = 0; i < 4; i++)
	{
		starvec[i] = 0.0L;
		for(int j = 0; j < 4; j++)
			starvec[i] += startmp[values[i][j]] * starcvec[j] * signs[i][j];
	}

	return EXIT_SUCCESS;
}

int star_polar2cartesian(double vec[4], double right_ascension, double declination, double r)
{
	vec[0] = r;
	vec[1] = r * cosl(right_ascension)*cosl(declination);
	vec[2] = r * sinl(right_ascension)*cosl(declination);
	vec[3] = r * sinl(declination);

	return EXIT_SUCCESS;
}

int star_cartesian2polar(double *right_ascension, double *declination, double vec[4], int equatorial)
{
	double r = 0.0L;
	r = sqrt(pow(vec[1], 2) + pow(vec[2],2) + pow(vec[3],2));
	*declination =  acosl(vec[3] / r);
	*right_ascension = atan2l(vec[2], vec[1]);
	if(equatorial)
		*declination = M_PI / 2.0L - *declination;

	return EXIT_SUCCESS;
}
