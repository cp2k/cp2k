/* A minimal wrapper for socket communication.

Copyright (C) 2013, Joshua More and Michele Ceriotti

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


Contains both the functions that transmit data to the socket and read the data
back out again once finished, and the function which opens the socket initially.
Can be linked to a FORTRAN code that does not support sockets natively.

Functions:
   error: Prints an error message and then exits.
   open_socket_: Opens a socket with the required host server, socket type and
      port number.
   write_buffer_: Writes a string to the socket.
   read_buffer_: Reads data from the socket.
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <sys/un.h>
#include <netdb.h>

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
      if (ai_err!=0) { perror("Error fetching host data. Wrong host name?"); exit(-1); }

      // creates socket
      sockfd = socket(res->ai_family, res->ai_socktype, res->ai_protocol);
      if (sockfd < 0) { perror("Error opening socket"); exit(-1); }
    
      // makes connection
      if (connect(sockfd, res->ai_addr, res->ai_addrlen) < 0) { perror("Error opening INET socket: wrong port or server unreachable"); exit(-1); }
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
      if (connect(sockfd, (struct sockaddr *) &serv_addr, sizeof(serv_addr)) < 0) { perror("Error opening UNIX socket: path unavailable, or already existing"); exit(-1); }
   }


   *psockfd=sockfd;
}

void writebuffer(int *psockfd, char *data, int* plen)
/* Writes to a socket.

Args:
   psockfd: The id of the socket that will be written to.
   data: The data to be written to the socket.
   plen: The length of the data in bytes.
*/

{
   int n;
   int sockfd=*psockfd;
   int len=*plen;

   n = write(sockfd,data,len);
   if (n < 0) { perror("Error writing to socket: server has quit or connection broke"); exit(-1); }
}


void readbuffer(int *psockfd, char *data, int* plen)
/* Reads from a socket.

Args:
   psockfd: The id of the socket that will be read from.
   data: The storage array for data read from the socket.
   plen: The length of the data in bytes.
*/

{
   int n, nr;
   int sockfd=*psockfd;
   int len=*plen;

   n = nr = read(sockfd,data,len);

   while (nr>0 && n<len )
   {  nr=read(sockfd,&data[n],len-n); n+=nr; }

   if (n == 0) { perror("Error reading from socket: server has quit or connection broke"); exit(-1); }
}


