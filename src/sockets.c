#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/stat.h>
#include <netinet/in.h>
#include <sys/un.h>
#include <netdb.h> 


//#define FS_WAIT 1     // uncomment to use a file system lock rather than mpi_barrier.

void error(const char *msg)
{   perror(msg);   }

void open_socket_(int *psockfd, int* inet, int* port, char* host)  // the darn fortran passes an extra argument for the string length. here I just ignore it
{
   int sockfd, portno, n;
   struct hostent *server;

   fprintf(stderr, "Connection requested %s, %d, %d\n", host, *port, *inet);
   struct sockaddr * psock; int ssock;
   if (*inet!=0)
   {  
      struct sockaddr_in serv_addr;      psock=(struct sockaddr *)&serv_addr;     ssock=sizeof(serv_addr);
      sockfd = socket(AF_INET, SOCK_STREAM, 0);
      if (sockfd < 0)  error("ERROR opening socket");
   
      server = gethostbyname(host);
      if (server == NULL)
      {
         fprintf(stderr, "ERROR, no such host %s \n", host);
         exit(-1);
      }

      bzero((char *) &serv_addr, sizeof(serv_addr));
      serv_addr.sin_family = AF_INET;
      bcopy((char *)server->h_addr, (char *)&serv_addr.sin_addr.s_addr, server->h_length);
      serv_addr.sin_port = htons(*port);
   }
   else
   {
      struct sockaddr_un serv_addr;      psock=(struct sockaddr *)&serv_addr;     ssock=sizeof(serv_addr);
      sockfd = socket(AF_UNIX, SOCK_STREAM, 0);
      bzero((char *) &serv_addr, sizeof(serv_addr));
      serv_addr.sun_family = AF_UNIX;
      strcpy(serv_addr.sun_path, "/tmp/ipi_");
      strcpy(serv_addr.sun_path+9, host);
   }
   
   if (connect(sockfd, psock, ssock) < 0) error("ERROR connecting");

   *psockfd=sockfd;
}

void writebuffer_(int *psockfd, char *data, int* plen)
{
   int n;   
   int sockfd=*psockfd;
   int len=*plen;

   n = write(sockfd,data,len);
   if (n < 0) error("ERROR writing to socket\n");
}


int readbuffer_(int *psockfd, char *data, int* plen)
{
   int n, nr;
   int sockfd=*psockfd;
   int len=*plen;

   n = nr = read(sockfd,data,len);
   
   while (nr>0 && n<len )
   {  nr=read(sockfd,&data[n],len-n); n+=nr; }
   if (n <= 0)  error("ERROR reading from socket\n"); 

   return n;
}

int check_reg(const char *path) {
       struct stat sb;
       return stat(path, &sb) == 0 && S_ISREG(sb.st_mode);
}

#ifdef FS_LOCK
int slock_(int *node, int *ionode)
{
//  fprintf(stderr,"creating file lock\n");
  if ((*node)==(*ionode))
  {    
    FILE *fh = fopen(".fs_sync", "w");
    fclose(fh);  
  }
}

int swait_(int *usec, int *node, int *ionode)
{
 //  fprintf(stderr, "swait %d %d %d %d\n", *node, *ionode, check_reg(".fs_sync"), ((*node)==(*ionode)) );
   if ((*node)==(*ionode))
   {  unlink(".fs_sync"); }
   else while(check_reg(".fs_sync")) {  usleep(*usec);    }
}

#else 
// just do nothing
int slock_(int *node, int *ionode) { return 0; }
int swait_(int *usec, int *node, int *ionode) { return 0; }
#endif

