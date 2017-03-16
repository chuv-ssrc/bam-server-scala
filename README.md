
![Travis build status](https://travis-ci.org/jdelafon/bam-server-scala.svg "Travis build status")


A simple BAM server in Scala
============================

The application described in this document is a microservice that allows querying portions of BAM files potentially located on a remote server, and implementing the OAuth2 authorization protocol. We name it "bam-server".

The goal is to provide

1. A way to query remote BAM files securely, according to user permissions;
2. An example widget that displays read alignments the user can access to;
3. Compatibility with any OAuth2-based authorization server.

Bam-server is written in Scala using the Play framework, and in itself represents only
the resource server in the OAuth2 diagram. It communicates with a local MySQL database
that maps each BAM file name with the corresponding sample identifier.

Functionalities
===============

Bam-server can return portions of BAM files using 3 different strategies:

1. Extracting a range of bytes;
2. Using samtools if available;
3. Using the Picard-tools library.

REST API
========

Let `:sample` be a sample identifier in the database;

Let `:token` be a Bearer token (JWT) for authentication;

Let `<range>` be a bytes range, formatted as "123-456";

Let `<region>` be a genomic region, formatted as "chr1:10000-20000".

    GET /

Says "BAM server operational.", just to test if the server is listening.

    POST /bai
    GET /bai/:sample/:token

Returns the content of the index (.bai) for that BAM file.

    POST /bam/range
    GET /bam/range/:sample/:token?range=<range>

Returns the content of the BAM file, expecting a Range HTTP header
to extract only the bytes range of interest - likely based on the index.

    POST /bam/samtools?region=<region>
    GET /bam/samtools/:sample/:token?region=<region>

Uses samtools (if available) to extract the region (``samtools view -hb <bam> <region>``).
Return the content as binary.

    POST /bam/json?region=<region>
    GET /bam/json/:sample/:token?region=<region>

Returns the reads for the given region in JSON format, using the [htsjdk](http://samtools.github.io/htsjdk/) library.
The fields correspond to the SAM file columns:

    [
      {
         "name": "HISEQ:206:C8E95ANXX:3:2113:2451:6639",   // read name
         "flag": 99,
         "chrom": "chr1",       // reference name
         "start": 1234,         // leftmost mapping position
         "end": 1334,           // rightmost mapping position
         "mapq": 50,            // mapping quality
         "cigar": "101M",       // cigar string
         "rnext": "=",          // 
         "pnext": 4567,         //
         "tlen": 283,           // template length, aka insert size
         "seq": "AATTAGGA...",  // [ACGTN=.],
         "qual": "AB<B@G>F..."  // per-base quality
      },
    ...
    ]

Configuration
=============

Settings can be edited in `conf/application.conf`. In particular,

- You must indicate ``env.BAM_PATH`` to tell where it can find BAM files for reading.
- It is set to use MySQL. You can change the database driver in the ``db`` section.
- You will probably want to configure CORS in ``play.filters.cors``.

Development
===========

Run the local dev server:

    activator run

Run tests:

    activator test


Deployment in production
========================

Build the source (to `target/universal/`):

    activator dist

Copy the built archive to destination, decompress.

Launch the server:

    ./bin/bam-server -v \
        -Dconfig.file=$SETTINGS \
        -Dhttp.port=$PORT \
        -Dhttps.port=$PORT_HTTPS

Test that it works:

    curl -k http://localhost:9000/

Configure a proxy to make it available from the outside world. Here with Apache:

    <VirtualHost *:443>
        ...
        ProxyPass         /bamserver  http://localhost:9000
        ProxyPassReverse  /bamserver  http://localhost:9000
        ...
    </Virtualhost>

The service becomes available at ``https://<host>/bamserver/``

N.B. To use another proxy than Apache, see
[Play HTTPServer docs](https://www.playframework.com/documentation/2.5.x/HTTPServer).

