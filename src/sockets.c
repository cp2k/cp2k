/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2021 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*  Copyright (C) 2013, Joshua More and Michele Ceriotti                      */
/*                                                                            */
/*  Permission is hereby granted, free of charge, to any person obtaining     */
/*  a copy of this software and associated documentation files (the           */
/*  "Software"), to deal in the Software without restriction, including       */
/*  without limitation the rights to use, copy, modify, merge, publish,       */
/*  distribute, sublicense, and/or sell copies of the Software, and to        */
/*  permit persons to whom the Software is furnished to do so, subject to     */
/*  the following conditions:                                                 */
/*                                                                            */
/*  The above copyright notice and this permission notice shall be included   */
/*  in all copies or substantial portions of the Software.                    */
/*                                                                            */
/*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,           */
/*  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF        */
/*  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.    */
/*  IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY      */
/*  CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,      */
/*  TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE         */
/*  SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                    */
/*----------------------------------------------------------------------------*/

/*******************************************************************************
 * \brief A minimal wrapper for socket communication.
 *        Contains both the functions that transmit data to the socket and read
 *        the data back out again once finished, and the function which opens
 *        the socket initially. Can be linked to a FORTRAN code that does not
 *        support sockets natively.
 * \author Joshua More and Michele Ceriotti
 ******************************************************************************/
#ifndef __NO_IPI_DRIVER

#define _POSIX_C_SOURCE 200809L

#include <math.h>
#include <netdb.h>
#include <netinet/in.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <sys/un.h>
#include <time.h>
#include <unistd.h>

/*******************************************************************************
 * \brief Opens a socket.
 * \param psockfd The id of the socket that will be created.
 * \param inet    An integer that determines whether the socket will be an inet
 *                or unix domain socket. Gives unix if 0, inet otherwise.
 * \param port    The port number for the socket to be created. Low numbers are
 *                often reserved for important channels, so use of numbers of 4
 *                or more digits is recommended.
 * \param host    The name of the host server.
 * \note  Fortran passes an extra argument for the string length, but this is
 *        ignored here for C compatibility.
 ******************************************************************************/
void open_socket(int *psockfd, int *inet, int *port, char *host) {
  int sockfd, ai_err;

  if (*inet > 0) { // creates an internet socket

    // fetches information on the host
    struct addrinfo hints, *res;
    char service[256];

    memset(&hints, 0, sizeof(hints));
    hints.ai_socktype = SOCK_STREAM;
    hints.ai_family = AF_INET;
    hints.ai_flags = AI_PASSIVE;

    sprintf(service, "%d", *port); // convert the port number to a string
    ai_err = getaddrinfo(host, service, &hints, &res);
    if (ai_err != 0) {
      perror("Error fetching host data. Wrong host name?");
      exit(-1);
    }

    // creates socket
    sockfd = socket(res->ai_family, res->ai_socktype, res->ai_protocol);
    if (sockfd < 0) {
      perror("Error opening socket");
      exit(-1);
    }

    // makes connection
    if (connect(sockfd, res->ai_addr, res->ai_addrlen) < 0) {
      perror("Error opening INET socket: wrong port or server unreachable");
      exit(-1);
    }
    freeaddrinfo(res);
  } else { // creates a unix socket
    struct sockaddr_un serv_addr;

    // fills up details of the socket address
    memset(&serv_addr, 0, sizeof(serv_addr));
    serv_addr.sun_family = AF_UNIX;
    strcpy(serv_addr.sun_path, "/tmp/ipi_");
    strcpy(serv_addr.sun_path + 9, host);

    // creates the socket
    sockfd = socket(AF_UNIX, SOCK_STREAM, 0);

    // connects
    if (connect(sockfd, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0) {
      perror(
          "Error opening UNIX socket: path unavailable, or already existing");
      exit(-1);
    }
  }

  *psockfd = sockfd;
}

/*******************************************************************************
 * \brief Writes to a socket.
 * \param psockfd The id of the socket that will be written to.
 * \param data    The data to be written to the socket.
 * \param plen    The length of the data in bytes.
 ******************************************************************************/
void writebuffer(int *psockfd, char *data, int *plen) {
  int n;
  int sockfd = *psockfd;
  int len = *plen;

  n = write(sockfd, data, len);
  if (n < 0) {
    perror("Error writing to socket: server has quit or connection broke");
    exit(-1);
  }
}

/*******************************************************************************
 * \brief Reads from a socket.
 * \param psockfd The id of the socket that will be read from.
 * \param data    The storage array for data read from the socket.
 * \param plen    The length of the data in bytes.
 ******************************************************************************/
void readbuffer(int *psockfd, char *data, int *plen) {
  int n, nr;
  int sockfd = *psockfd;
  int len = *plen;

  n = nr = read(sockfd, data, len);

  while (nr > 0 && n < len) {
    nr = read(sockfd, &data[n], len - n);
    n += nr;
  }

  if (n == 0) {
    perror("Error reading from socket: server has quit or connection broke");
    exit(-1);
  }
}

/*******************************************************************************
 * \brief Mini-wrapper to nanosleep
 * \param dsec number of seconds to wait (float values accepted)
 ******************************************************************************/
void uwait(double *dsec) {
  struct timespec wt, rem;
  wt.tv_sec = floor(*dsec);
  wt.tv_nsec = (*dsec - wt.tv_sec) * 1000000000;
  nanosleep(&wt, &rem);
}

#endif
