
A simple BAM server in Scala
============================

The application described in this document is a microservice that allows querying portions of BAM files potentially located on a remote server, and implementing the OAuth2 authorization protocol. We name it "bam-server".

The goal is to provide

1. A way to query remote BAM files securely, according to user permissions;
2. An example widget that displays read alignments the user can access to;
3. Compatibility with any OAuth2, token-based authorization server.

Bam-server is written in Scala using the Play framework, and in itself represents only
the resource server in the OAuth2 diagram. It communicates with a local database
that maps each BAM file name with the corresponding sample identifier.

Bam-server can return portions of BAM files using 3 different strategies:

1. Extracting a range of bytes;
2. Using samtools if available;
3. Using the Picard-tools library.

For the complete documentation, see

[](http://bam-server-scala.readthedocs.io/en/latest/src/about.html)


