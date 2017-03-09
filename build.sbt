import NativePackagerHelper._

name := """bam-server"""
version := "1.0-SNAPSHOT"
scalaVersion := "2.11.7"

lazy val root = (project in file(".")).enablePlugins(PlayScala)

libraryDependencies ++= Seq(
  jdbc,
  cache,
  ws,
  filters,
  evolutions,
  "com.github.samtools" % "htsjdk" % "2.1.1"
)

// Testing
libraryDependencies ++= Seq(
  "org.scalatestplus.play" %% "scalatestplus-play" % "1.5.0" % Test,
  "org.mockito" % "mockito-all" % "1.9.5"
)

// Database drivers
libraryDependencies ++= Seq(
  "mysql" % "mysql-connector-java" % "5.1.36",
  "org.xerial" % "sqlite-jdbc" % "3.16.1",
  "com.h2database" % "h2" % "1.4.193"
)

// JWT
libraryDependencies ++= Seq(
  "com.pauldijou" %% "jwt-core" % "0.8.1",
  "com.pauldijou" %% "jwt-play" % "0.8.1",
  "org.mindrot" % "jbcrypt" % "0.3m"
)

resolvers += "scalaz-bintray" at "http://dl.bintray.com/scalaz/releases"

// Add non-default directories to build (/target/universal/) using NativePackageHelper
mappings in Universal ++= directory("scripts")

// Use application.test.conf for testing
javaOptions in Test += "-Dconfig.file=conf/application.test.conf"
