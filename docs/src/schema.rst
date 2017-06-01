

Local database
==============

Here we describe the local bam-server database.
It links together app, user and sample identifers to set individual permissions.


Schema
------

.. table::

    =============  ===========
    Table name     Description
    =============  ===========
    apps           The client applications, that query the bam-server.
    users          The individual users of an app, as identified by auth tokens coming from the app.
    samples        A sample corresponds to one BAM file.
    users_samples  The attribution of samples to users
    =============  ===========

Here is the complete schema creation script (for MySQL):

.. code:: MySQL

    CREATE TABLE `apps` (
        `id` INTEGER(11) NOT NULL PRIMARY KEY AUTO_INCREMENT,
        `iss` VARCHAR(255) NOT NULL,
        `key` TEXT NOT NULL,
        `algorithm` VARCHAR(255) DEFAULT 'RS256',
        `description` VARCHAR(255) DEFAULT NULL,
        `isActive` TINYINT(1) NOT NULL DEFAULT 1
    );

    CREATE TABLE `users` (
        `id` INTEGER(11) NOT NULL PRIMARY KEY AUTO_INCREMENT,
        `app_id` INTEGER(11) NOT NULL,
        `username` VARCHAR(255) NOT NULL,
        `group` VARCHAR(255) DEFAULT NULL,
        `isActive` TINYINT(1) DEFAULT 1,
        `isAdmin` TINYINT(1) DEFAULT 0,
        FOREIGN KEY (`app_id`) REFERENCES `apps`(`id`)
    );

    CREATE TABLE `samples` (
        `id` INTEGER(11) NOT NULL PRIMARY KEY AUTO_INCREMENT,
        `name` VARCHAR(255) NOT NULL,
        `filename` VARCHAR(255) NOT NULL,
        `project` VARCHAR(255) DEFAULT NULL,
        `hash` VARCHAR(255) DEFAULT NULL,
        `description` VARCHAR(255) DEFAULT NULL,
        `isOndisk` TINYINT(1) DEFAULT NULL,
        `isActive` TINYINT(1) NOT NULL DEFAULT 1
    );

    CREATE TABLE `users_samples` (
        `id` INTEGER(11) NOT NULL PRIMARY KEY AUTO_INCREMENT,
        `user_id` INTEGER(11) NOT NULL,
        `sample_id` INTEGER(11) NOT NULL,
        `isActive` TINYINT(1) NOT NULL DEFAULT 1,
        FOREIGN KEY (`user_id`) REFERENCES `users`(`id`),
        FOREIGN KEY (`sample_id`) REFERENCES `samples`(`id`)
    );


Fields description
------------------

.. table::

    =========================== ===========
    **apps.iss**                App identifier - the "iss" claim of JWTs.
    **apps.key**                The token signature verification key (a public key for RSA, a shared secret for HMAC).
    **apps.algorithm**          The signature algorithm (one of HS256,HS384,HS512,RS256,RS384,RS512).
    **apps.description**        [optional] Some text to describe the app.
    **users.app_id**            The app ID.
    **users.username**          The user identifier - the one passed in the JWT to identify the requester.
    **users.group**             [optional] The name of a group the user belong to :sup:`1`.
    **users.isAdmin**           Whether the user has administrator rights (can create or delete users, samples, etc.).
    **samples.name**            The sample identifier - the one used in the API to query the corresponding BAM file.
    **samples.filename**        The name of the corresponding BAM file in ``env.BAM_PATH`` (see Configuration).
    **samples.project**         [optional] The name of a project the sample belongs to :sup:`1`.
    **samples.hash**            [optional] The hash of the BAM file :sup:`1`.
    **samples.description**     [optional] Some text to describe the sample.
    **samples.isOnDisk**        [optional] To mark a file as not found on disk anymore :sup:`1`.
    **users_samples.user_id**   The user ID.
    **users_samples.sample_id** The sample ID.
    =========================== ===========

[1] These fields are not used directly, but can be useful for automatic management tasks.


Some test data
--------------

For convenience, this generates 3 apps with different signature algorithms,
2 users (1 admin), 2 samples, and attributes samples to the users.

.. code:: MySQL

    INSERT INTO `apps` VALUES
    (1, 'https://jdelafon.eu.auth0.com/', '-----BEGIN RSA PUBLIC KEY-----\nMIIBIjANBgkqhkiG9w0BAQEFAAOCAQ8AMIIBCgKCAQEAtWWKxPv9vsWdRR/hcmJF\nsQjjUrMs/OsVstyNJXwmWuhl3lNIZwwEDoJbnE9IKPyizyNwbnB9FmJnClCboUeP\nbkuIrDM63+S+PtX/SQ9YI5yDxz+88dRYT86WP23wcWMO3txV2GAu62RVGSl48ZJP\nSyu94NBIiZOO5oDJpWDInhZphiMQ3u/rEwlVxVMt0CTTInfl4iX0sCtymD2y6M38\nVrQwHOzSddFrbI58t4Rfal4SttwdmXONRnj7mrgl5G6v7IHEa/HOrlT1rSLOMBKz\nOfmZy+bdlt5zrx3Adfzgn1BC6DGlG3Y9QYMOPpXjbzRO3rv9Fl5bRJyn5Ih82Cey\ndQIDAQAB\n-----END RSA PUBLIC KEY-----', 'RS256', 'Using a public certificate in .cer format', 1),
    (2, 'testapp', '-----BEGIN RSA PUBLIC KEY-----\nMFwwDQYJKoZIhvcNAQEBBQADSwAwSAJBAK7ttYaE/1ldsb0OJQDQhhDWqwuFWIyt\nxgYIJH1HYA4UpA/Nm24fERIA1xi2Pomep6VTnQ/ThFP5hn2NyITwCIsCAwEAAQ==\n-----END RSA PUBLIC KEY-----', 'RS256', 'Using a public key in .pem format', 1),
    (3, 'testapp-hmac', 'secretHMACkey', 'HS256', 'Using a shared secret key', 1);

    INSERT INTO `users`(`id`,`app_id`,`username`,`isActive`,`isAdmin`) VALUES
    (1, 1, 'admin@test.com', 1, 1),
    (2, 1, 'test@test.com', 1, 0),
    (3, 2, 'test@test.com', 1, 0),
    (4, 3, 'test@test.com', 1, 0);

    INSERT INTO `samples`(`id`,`name`,`filename`,`isActive`) VALUES
    (1, 'sample1', 'test1.bam', 1),
    (2, 'sample2', 'test2.bam', 1);

    INSERT INTO `users_samples`(`user_id`,`sample_id`) VALUES
    (1, 1),(1, 2),(2, 1);

