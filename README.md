
![Travis build status](https://travis-ci.org/jdelafon/bam-server-scala.svg "Travis build status")


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

Functionalities
===============

Bam-server can return portions of BAM files using 3 different strategies:

1. Extracting a range of bytes;
2. Using samtools if available;
3. Using the Picard-tools library.

REST API
========

Let `:sample`/`<sample>` be a sample identifier in the database;

Let `<range>` be a bytes range, formatted as "123-456";

Let `<region>` be a genomic region, formatted as "chr1:10000-20000".

    GET /

Says "BAM server operational.", just to test if the server is listening.

    POST /bai                                   # data = '{"sample": "<sample>"}'
    GET /bai/:sample

Returns the content of the index (.bai) for that BAM file.

    POST /bam/range                             # data = '{"sample": "<sample>"}'
    GET /bam/range/:sample

Returns the content of the BAM file, expecting a Range HTTP header
to extract only the bytes range of interest - likely based on the index.
The bytes range can also be passed as an argument to the request (`?range=<range>`).
Return the content as binary.

    POST /bam/samtools?region=<region>          # data = '{"sample": "<sample>"}'
    GET /bam/samtools/:sample?region=<region>

Uses samtools (if available) to extract the region (``samtools view -hb <bam> <region>``).
Return the content as binary.

    POST /bam/json?region=<region>              # data = '{"sample": "<sample>"}'
    GET /bam/json/:sample?region=<region>

Returns the reads for the given region in JSON format, 
using the [htsjdk](http://samtools.github.io/htsjdk/) library.
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
    
All requests expect an Authorization header to be added to the request,
with a RSA256-signed Bearer token. For example:

    curl -i -H "Authorization: Bearer xxxx.yyyy.zzzz" -X POST -d '{"sample": "SAMPLE1"} http://localhost:9000/bai
    curl -i -H "Authorization: Bearer xxxx.yyyy.zzzz" http://localhost:9000/bai/SAMPLE1

Alternatively, all request can accept a URL argument `?token=<JWT>` to pass the token.
This mode is not recommended but sometimes necessary, depending on client-side constraints.
For example:

    curl -i http://localhost:9000/bai/SAMPLE1?token=xxxx.yyyy.zzzz

For more details, see Authorization below.

Authorization
=============

To ensure that only registered users have access only to the data that belongs to them,
all requests must include a Json Web Token (JWT) in their "Authorization" header.
It is a **RSA256-signed** Bearer token which payload includes at least the user
identifer under claim `"name"` and the app identifier under claim `"iss"`. For example

    # Header
    {
      "typ": "JWT",
      "alg": "RS256"
    }
    
    # Payload 
    {
      "name": "myUsername77",
      "iss": "myAppname",
      "exp": 31490863741,
      "iat": 1490863741,
      "sub": ...,
      "aud": ...,
      ...
    }
    
    # Signature
    RSASHA256(
      base64UrlEncode(header) + "." +
      base64UrlEncode(payload),
      <public key>,
      <private key>
    )
    
    # Final token
    <Base64(Header)>.<Base64(Payload)>.<Base64(Signature)>
    
For more details on JWTs, see [jwt.io](jwt.io).

The client app/authorization server encrypts the JWT with its secret RSA private key.
To verify it, the public key or public certificate (.cer, .pem) must be copied into 
`resources/rsa_keys`. The `apps` table matches the "iss" claim with the name of the
public certificate.

Public keys in PEM format look like this:

> -----BEGIN RSA PUBLIC KEY-----
> MFwwDQYJKoZIhvcNAQEBBQADSwAwSAJBAK7ttYaE/1ldsb0OJQDQhhDWqwuFWIyt
> xgYIJH1HYA4UpA/Nm24fERIA1xi2Pomep6VTnQ/ThFP5hn2NyITwCIsCAwEAAQ==
> -----END RSA PUBLIC KEY-----

Certificates in CER format look like this:

> -----BEGIN CERTIFICATE-----
> MIIC8jCCAdqgAwIBAgIJY0fhkAkJqc1YMA0GCSqGSIb3DQEBBQUAMCAxHjAcBgNV
> ...
> -----END CERTIFICATE-----

The authorization process then works as follows:

1. The user ("name") is matched against the `users` table in the database
   to make sure he is registered.

2. The app name ("iss") is matched against the `apps` table in the database,
   so as to retreive the corresponding RSA public key file. 
   The `resources/rsa_keys/` directory is scanned recursively for .cer/.pem files,
   and when one is found with the right name, it is parsed to get the public key.
  
3. The JWT is verified using the public key.

Configuration
=============

Settings can be edited in `/conf/application.conf`. In particular,

- You must indicate ``env.BAM_PATH`` to tell where it can find BAM files for reading.
- It is set to use SQLite. You can change the database driver in the ``db`` section.
- You will probably want to configure CORS in ``play.filters.cors``.

A default SQLite database is provided with the project source ("/resources/db/bam-server")

Development
===========

Run the local dev server:

    activator run

Run tests:

    activator test
    
Tests run on a temporary, in-memory H2 database, 
by running evolutions (/conf/evolutions/default/1.sql) from a clean state.


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

Configure a proxy to make it available from the outside world. Example with Apache:

    <VirtualHost *:443>
        ...
        ProxyPass         /bamserver  http://localhost:9000
        ProxyPassReverse  /bamserver  http://localhost:9000
        ...
    </Virtualhost>

The service becomes available at ``https://<host>/bamserver/``

N.B. To use another proxy than Apache, see
[Play HTTPServer docs](https://www.playframework.com/documentation/2.5.x/HTTPServer).

