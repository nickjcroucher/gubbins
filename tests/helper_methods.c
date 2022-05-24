#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include "helper_methods.h"


int compare_files(char expected_output_filename[],char actual_output_filename[] )
{
  FILE *expected_output_fh;
  FILE *actual_output_fh;
  
  char    *expected_buffer;
  char    *actual_buffer;
  long    numbytes;
  
  expected_output_fh = fopen(expected_output_filename, "r");
  actual_output_fh = fopen(actual_output_filename, "r");
  
  fseek(expected_output_fh, 0L, SEEK_END);
  numbytes = ftell(expected_output_fh);
  fseek(expected_output_fh, 0L, SEEK_SET);	
  expected_buffer = (char*)calloc(numbytes +1, sizeof(char));	
  fread(expected_buffer, sizeof(char), numbytes, expected_output_fh);
  fclose(expected_output_fh);
  
  fseek(actual_output_fh, 0L, SEEK_END);
  numbytes = ftell(actual_output_fh);
  fseek(actual_output_fh, 0L, SEEK_SET);	
  actual_buffer = (char*)calloc(numbytes +1, sizeof(char));	
  fread(actual_buffer, sizeof(char), numbytes, actual_output_fh);
  fclose(actual_output_fh);
  
  if(strcmp(expected_buffer,actual_buffer) == 0)
  { 
    free(expected_buffer);
    free(actual_buffer);
    return 1;
  }

  free(expected_buffer);
  free(actual_buffer);
  
  return 0;
}

int number_of_recombinations_in_file(char * fileName)
{
	FILE *file = fopen(fileName, "r");
	int ch, prev = '\n', lines = 0;
	while ( (ch = fgetc(file)) != EOF )
	{
		if ( ch == '\n' )
		{
			++lines; 
		}
		prev = ch; 
	}
	fclose(file);
	if ( prev != '\n' ) 
	{
		++lines;
	}
  return  (int) lines/6;
}

int file_exists(char * fileName)
{
   struct stat buf;
   int i = stat ( fileName, &buf );
     if ( i == 0 )
     {
       return 1;
     }
     return 0;
}



int cp(const char *to, const char *from)
{
    int fd_to, fd_from;
    char buf[4096];
    ssize_t nread;
    int saved_errno;

    fd_from = open(from, O_RDONLY);
    if (fd_from < 0)
        return -1;

    fd_to = open(to, O_WRONLY | O_CREAT | O_EXCL, 0666);
    if (fd_to < 0)
        goto out_error;

    while (nread = read(fd_from, buf, sizeof buf), nread > 0)
    {
        char *out_ptr = buf;
        ssize_t nwritten;

        do {
            nwritten = write(fd_to, out_ptr, nread);

            if (nwritten >= 0)
            {
                nread -= nwritten;
                out_ptr += nwritten;
            }
            else if (errno != EINTR)
            {
                goto out_error;
            }
        } while (nread > 0);
    }

    if (nread == 0)
    {
        if (close(fd_to) < 0)
        {
            fd_to = -1;
            goto out_error;
        }
        close(fd_from);

        /* Success! */
        return 0;
    }

  out_error:
    saved_errno = errno;

    close(fd_from);
    if (fd_to >= 0)
        close(fd_to);

    errno = saved_errno;
    return -1;
}
