
A simple BAM server in Scala
============================

Bam-server is a microservice that allows querying slices of BAM files remotely and securely.

Bam-server is written in Scala using the Play framework, and in itself represents only
the resource server in the OAuth2 diagram. It communicates with a local database
that maps each BAM file name with the corresponding sample identifier.

Bam-server can return portions of BAM files using 3 different strategies:

1. Extracting a range of bytes;
2. Using samtools if available to return reads form a region in BAM format;
3. Using htsjdk to return reads as JSON.

For the complete documentation, see

[http://bam-server-scala.readthedocs.io/en/latest/src/about.html](http://bam-server-scala.readthedocs.io/en/latest/src/about.html)


