
![Travis build status](https://travis-ci.org/jdelafon/bam-server-scala.svg "Travis build status")


A simple BAM server in Scala
============================

The application described in this document is a microservice that allows querying portions of BAM files potentially located on a remote server, and implementing the OAuth2 authorization protocol. We name it "bam-server".

The goal is to provide

1. A way to query remote BAM files securely, according to user permissions;
2. An example widget that displays read alignments the user can access to;
3. Compatibility with any OAuth2-based authorization server.

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

Let `:sample` be a sample identifier in the database;

Let `<range>` be a bytes range, formatted as "123-456";

Let `<region>` be a genomic region, formatted as "chr1:10000-20000".

    GET /

Says "BAM server operational.", just to test if the server is listening.

    POST /bai
    GET /bai/:sample

Returns the content of the index (.bai) for that BAM file.

    POST /bam/range
    GET /bam/range/:sample?range=<range>

Returns the content of the BAM file, expecting a Range HTTP header
to extract only the bytes range of interest - likely based on the index.

    POST /bam/samtools?region=<region>
    GET /bam/samtools/:sample?region=<region>

Uses samtools (if available) to extract the region (``samtools view -hb <bam> <region>``).
Return the content as binary.

    POST /bam/json?region=<region>
    GET /bam/json/:sample?region=<region>

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


Example of usage with Auth0
===========================

[Auth0](https://manage.auth0.com) is a commercial solution that provides an easy-to-manage authorization server
(for free if you have few users and usual requirements). 
Its only role is to identify users of the client app and return valid tokens. 

1. First time setup 

    1. Create an Auth0 account

    2. Create a new client app, say "bam-server". It is given a Client ID and a Client secret.
    
    3. Go to "Clients > Show advanced settings". In the "OAuth" tab, choose "RS256". This changes
       the encryption protocol from the default HMAC (symmetric, shared secret) 
       to RSA (asymmetric private/public - more secure).
       
    4. Still in the advanced settings, copy the public certificate from "Certificates > Signing Certificate",
       or download it in CER/PEM format.
    
    5. Copy the certificate to "test.cer" in the `resources/rsa_keys` directory. 
       So now you have a file `resources/rsa_keys/test.cer`. 

2. Obtain a token (for testing; the client app is out of this scope)

    1. Let's say your client web app is served as localhost:3000.
    
    2. Add "localhost:3000" to Auth0 "Clients > Allowed callback URLs".
    
    3. Point your browser to 
         
           https://<domain>/authorize
               ?scope=openid name nickname email user_id
               &response_type=token
               &client_id=<client_id>
               &redirect_uri=localhost:3000
               &nonce=1234&state=1234

       where you replace `<domain>` and `<client_id>` by your own Auth0 client settings.
       N.B. The "scope" argument is what permits to have the "name" field in the token.
       
    4. Log in; you get redirected to localhost:3000/ and you can find an "id_token" in the url.
       This is a JWT that bam-server accepts (*not* "access_token"). For instance:
       
       > http://<i></i>localhost:3000/#access_token=ZVvVzzBIucYdhnW3&expires_in=86400&id_token=**xxxxxxx.yyyyyyy.zzzzzzz**&token_type=Bearer&state=1234  

3. Prepare the database

    1. Tell your application where to find the certificate by adding a row to the `apps` table
       of the database:
       
           INSERT INTO apps(iss, name) VALUES ('https://<auth0_domain>/', 'test', NULL);
           
       where you replace `<auth0_domain>` by your own Auth0 Domain setting, and `'test'` refers to the 
       public key file you just created.
       N.B. A typical "iss" claim for Auth0 looks like "https://<domain>/", but in principle it can be anything else.
       
       N.B. This functionality will become part of the REST API.
       
    2. Add a user for this app to the database:
    
           INSERT INTO users(app_id, username) VALUES (<app_id>, <username>);
           
       where the `<app_id>` is the id of the row inserted above, and `<username>` corresponds
       to the "name" claim of the JWT (the user name in the client).
       In Auth0, the `<username>` is your email.

       N.B. This functionality will become part of the REST API.
       
    3. Add a sample:
     
           INSERT INTO bam(sample, filename) VALUES (<sample_name>, <bam_filename>);
           
       where you replace `<sample_name>` by the unique identifier of the sample,
       and `<bam_filename>` by the name of the BAM file.
       
       Make sure the BAM file is available at `env.BAM_PATH` (see "Configuration" above).
    
       N.B. This functionality will become part of the REST API.
       
    4. Attribute the sample to a user:
    
           INSERT INTO users_bam(user_id, bam_id) VALUES (<user_id>, <bam_id>);
           
       where `<user_id>` and `<bam_id>` are the respective ids of the rows inserted above.

       N.B. This functionality will become part of the REST API.

3. Use the token

   Try to read a BAM index with curl:
    
       curl -i -H "Authorization: Bearer <token>" http://localhost:9000/bai/<sample_name>