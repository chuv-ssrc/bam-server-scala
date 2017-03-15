
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

Let `$key` be a sample identifier in the database.
Let `token` be a Bearer token (JWT) for authentication.


    GET /

Says "BAM server operational.", just to test if the server is listening.


    POST /bai
    GET /bai/:token

Returns the content of the index (.bai) for that BAM file.


    POST /bam/range
    GET /bam/range/:range/:token

Returns the content of the BAM file, expecting a Range HTTP header
to extract only the bytes range of interest - likely based on the index.


    POST /bam/samtools/:region
    GET /bam/samtools/:region/:token

Uses samtools (if available) to extract the region (``samtools view -hb <bam> <region>``).
Return the content as binary.


    POST /bam/json/
    GET /bam/json/:region/:token

Returns the reads for the given region in JSON format.
Currently, it looks like this:

    [
      {
         "chrom": "chr1",
         "name": "long-read-name",
         "start": 1234,
         "end": 45678,
         "cigar": "long-cigar-string"
      },
    ...
    ]

It can be edited to return whatever you want for you own BAM viewer.
Look at [htsjdk docs under SAMRecord](https://samtools.github.io/htsjdk/javadoc/htsjdk)
for available fields.

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


Deployment
==========

Build the source (to `target/universal/`):

    activator dist

Copy the built archive to destination, decompress.

Launch the server (see also `/scripts/deploy.sh`):

    ./bin/bam-server -v \
        -Dconfig.file=$SETTINGS \
        -Dhttp.port=$PORT \
        -Dhttps.port=$PORT_HTTPS

Test that it works (otherwise file an issue!):

    curl -k http://localhost:9000/

Configure a proxy. Here with Apache:

    <VirtualHost *:443>
        ...
        ProxyPass         /bamserver  http://localhost:9000
        ProxyPassReverse  /bamserver  http://localhost:9000
        ...
    </Virtualhost>

The service will then get available at ``https://<host>/bamserver/``

N.B. To use another proxy than Apache, see
[Play HTTPServer docs](https://www.playframework.com/documentation/2.5.x/HTTPServer).

