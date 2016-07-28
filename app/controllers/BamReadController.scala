package controllers

import java.nio.file.{Path, Paths}
import javax.inject.Inject

import play.api.db.Database
import play.api.mvc.Action
import utils.Constants._

import sys.process._

/**
  * Return the BAM text as output by "samtools view -hb <filename> <region>".
  * It is encoded in base64 because transfers have to be ascii (?).
  * Subject to "java.lang.OutOfMemoryError: Java heap space".
  */
class BamReadController @Inject()(db: Database) extends BamController(db) {

  def read(key: String, region: String, description: Option[String]) = Action {
    val filename = getBamName(key)
    val bamPath: Path = Paths.get(BAM_PATH, filename)
    val command = s"samtools view -hb $bamPath $region" #| "base64"
    val res: String = command.!!
    Ok(res).as(HTML)

    //import play.api.libs.streams.Streams
    //import akka.util.ByteString
    //import java.io.ByteArrayInputStream

    //val bam = bamPath.toFile
    //val stream: FileInputStream = new FileInputStream(bam)
    //val source = StreamConverters.fromInputStream(() => stream)
    //val contentLength = bam.length
    //Result(
    //  header = ResponseHeader(OK, Map(
    //    CONTENT_TYPE -> BINARY,
    //    CONTENT_LENGTH -> contentLength.toString,
    //    CONNECTION -> "keep-alive"
    //  )),
    //  body = HttpEntity.Streamed(source, Some(contentLength), Some(BINARY))
    //)
  }
}

