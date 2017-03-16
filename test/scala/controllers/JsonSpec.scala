package scala.controllers

import org.scalatest.BeforeAndAfter
import org.scalatestplus.play._
import play.api.libs.json._
import play.api.test.Helpers._
import play.api.test._


/**
  * Test JsonController.
  * The test bam has coordinates in range chr1:761997-762551.
  */
class JsonSpec extends PlaySpec with OneAppPerSuite with BeforeAndAfter {

  val token = "asdf"
  val testkey = "testkey"
  val body: JsValue = Json.parse(s"""{"key": "$testkey"}""")
  val headers = (AUTHORIZATION -> s"Bearer $token")

  "JsonController" should {

    "provide reads in JSON format if a region is given (POST)" in {
      val request = FakeRequest(POST, "/bam/json?region=chr1:761997-762551").withJsonBody(body).withHeaders(headers)
      val response = route(app, request).get
      status(response) mustBe OK
      contentType(response) must be(Some(JSON))

      val content = contentAsJson(response)
      //println(content.head)
      val read1 = content.head
      (read1 \ "name").as[String] must equal("HISEQ:206:C8E95ANXX:3:2113:2451:6639")
      (read1 \ "flag").as[Int] must equal(99)
      (read1 \ "chrom").as[String] must equal("chr1")
      (read1 \ "start").as[Int] must equal(761997)
      (read1 \ "end").as[Int] must equal(762097)  // + 100
      (read1 \ "mapq").as[Int] must equal(50)
      (read1 \ "cigar").as[String] must equal("101M")
      //(read1 \ "rnext").as[String] must equal("=")
      //(read1 \ "pnext").as[Int] must equal(762179)
      (read1 \ "tlen").as[Int] must equal(283)
      (read1 \ "seq").as[String] must equal("CTACTGACGGTCAAGGCCTCCTCATTGTATTCTGTCCTCCATATCTCTGCTGATTCCCATTTTGTCTATTTCCATTTACCCCACTACTGCTTGCTCAGGTC")
      (read1 \ "qual").as[String] must equal("AB<B@G>FAF=E@BHFFFAEFAF@?>G><?=FAG=EFAEF@><>EAEAGFAG>>=EFF@>===G=EA<>==EF>>==<FFCF@FA;FAAFA=GFAD?B6;C")
    }

    "provide reads in JSON format if a region is given (GET)" in {
      val request = FakeRequest(GET, s"/bam/json/$testkey/$token?region=chr1:761997-762551")
      val response = route(app, request).get
      status(response) mustBe OK
      contentType(response) must be(Some(JSON))
    }

  }

}