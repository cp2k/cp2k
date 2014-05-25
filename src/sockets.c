#ifndef __NO_IPI_DRIVER
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <sys/un.h>
#include <netdb.h> 

void error(const char *msg)
// Prints an error message and then exits.
{   perror(msg);  exit(-1);   }

void open_socket(int *psockfd, int* inet, int* port, char* host)
/* Opens a socket.

Note that fortran passes an extra argument for the string length, but this is
ignored here for C compatibility.

Args:
   psockfd: The id of the socket that will be created.
   inet: An integer that determines whether the socket will be an inet or unix
      domain socket. Gives unix if 0, inet otherwise.
   port: The port number for the socket to be created. Low numbers are often
      reserved for important channels, so use of numbers of 4 or more digits is
      recommended.
   host: The name of the host server.
*/

{
   int sockfd, ai_err;
   
   if (*inet>0)
   {  // creates an internet socket

      // fetches information on the host      
      struct addrinfo hints, *res;  
      char service[256];
   
      memset(&hints, 0, sizeof(hints));
      hints.ai_socktype = SOCK_STREAM;
      hints.ai_family = AF_UNSPEC;
      hints.ai_flags = AI_PASSIVE;

      sprintf(service,"%d",*port); // convert the port number to a string
      ai_err = getaddrinfo(host, service, &hints, &res); 
      if (ai_err!=0) error("Error fetching host data. Wrong host name?");

      // creates socket
      sockfd = socket(res->ai_family, res->ai_socktype, res->ai_protocol);
      if (sockfd < 0)  error("Error opening socket");
    
      // makes connection
      if (connect(sockfd, res->ai_addr, res->ai_addrlen) < 0) error("Error opening INET socket: wrong port or server unreachable");
      freeaddrinfo(res);
   }
   else
   {  // creates a unix socket
      struct sockaddr_un serv_addr;    

      // fills up details of the socket addres
      memset(&serv_addr, 0, sizeof(serv_addr));
      serv_addr.sun_family = AF_UNIX;
      strcpy(serv_addr.sun_path, "/tmp/ipi_");
      strcpy(serv_addr.sun_path+9, host);
  
      // creates the socket
      sockfd = socket(AF_UNIX, SOCK_STREAM, 0);

      // connects
      if (connect(sockfd, (struct sockaddr *) &serv_addr, sizeof(serv_addr)) < 0) error("Error opening UNIX socket: path unavailable, or already existing");
   }

   *psockfd=sockfd;
}

void writebuffer(int *psockfd, char *data, int* plen)
{
   int n;   
   int sockfd=*psockfd;
   int len=*plen;

   n = write(sockfd,data,len);
   if (n < 0) error("Error writing to socket: server has quit or connection broke");
}

void writebuffer_i(int *psockfd, int *data, int* plen)
{ int llen=(*plen)*4; writebuffer(psockfd, (char*) data, &llen); }

void writebuffer_d(int *psockfd, double *data, int* plen)
{ int llen=(*plen)*8; writebuffer(psockfd, (char*) data, &llen); }

int readbuffer(int *psockfd, char *data, int* plen)
{
   int n, nr;
   int sockfd=*psockfd;
   int len=*plen;

   n = nr = read(sockfd,data,len);

   while (nr>0 && n<len )
   {  nr=read(sockfd,&data[n],len-n); n+=nr; }

   if (n <= 0) error("Error reading from socket: server has quit or connection broke");

   return n;
}

int readbuffer_i(int *psockfd, int *data, int* plen)
{ int llen=(*plen)*4; return readbuffer(psockfd, (char*) data, &llen); }

int readbuffer_d(int *psockfd, double *data, int* plen)
{ int llen=(*plen)*8; return readbuffer(psockfd, (char*) data, &llen); }

#endif
