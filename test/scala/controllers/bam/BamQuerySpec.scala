package scala.controllers.bam

import org.scalatestplus.play._
import play.api.libs.json._
import play.api.test.Helpers._

import scala.setup.WithToken


/**
  * Test common authorization control methods of users on samples, in the /bai case.
  */
class BamQuerySpec extends PlaySpec with OneAppPerSuite with WithToken {

  def testGetAccess(root: String, userType: String, sampleName: String, args: String = ""): Int = {
    val response = if (userType == "admin") {
      route(app, FakeAdminRequest(GET, s"/$root/$sampleName$args")).get
    } else {
      route(app, FakeAuthorizedRequest(GET, s"/$root/$sampleName$args")).get
    }
    status(response)
  }

  def testPostAccess(root: String, userType: String, sampleName: String, args: String = ""): Int = {
    val body: JsValue = Json.parse(s"""{"sample": "$sampleName"}""")
    val response = if (userType == "admin") {
      route(app, FakeAdminRequest(POST, s"/$root").withJsonBody(body)).get
    } else {
      route(app, FakeAuthorizedRequest(POST, s"/$root").withJsonBody(body)).get
    }
    status(response)
  }

  "BAM query controllers" should {

    "return a 200 if the user can access the sample, 403 otherwise (GET)" in {
      testGetAccess("bai", "", testSample1) mustBe OK
      testGetAccess("bai", "", testSample2) mustBe FORBIDDEN
      testGetAccess("bai", "admin", testSample1) mustBe OK
      testGetAccess("bai", "admin", testSample2) mustBe OK
      testGetAccess("bam/json", "", testSample1) mustBe OK
      testGetAccess("bam/json", "", testSample2) mustBe FORBIDDEN
      testGetAccess("bam/json", "admin", testSample1) mustBe OK
      testGetAccess("bam/json", "admin", testSample2) mustBe OK
      testGetAccess("bam/range", "", testSample1) mustBe OK
      testGetAccess("bam/range", "", testSample2) mustBe FORBIDDEN
      testGetAccess("bam/range", "admin", testSample1) mustBe OK
      testGetAccess("bam/range", "admin", testSample2) mustBe OK
      testGetAccess("bam/samtools", "", testSample1) mustBe OK
      testGetAccess("bam/samtools", "", testSample2) mustBe FORBIDDEN
      testGetAccess("bam/samtools", "admin", testSample1) mustBe OK
      testGetAccess("bam/samtools", "admin", testSample2) mustBe OK
    }

    "return a 200 if the user can access the sample, 403 otherwise (POST)" in {
      testPostAccess("bai", "", testSample1) mustBe OK
      testPostAccess("bai", "", testSample2) mustBe FORBIDDEN
      testPostAccess("bai", "admin", testSample1) mustBe OK
      testPostAccess("bai", "admin", testSample2) mustBe OK
      testPostAccess("bam/json", "", testSample1) mustBe OK
      testPostAccess("bam/json", "", testSample2) mustBe FORBIDDEN
      testPostAccess("bam/json", "admin", testSample1) mustBe OK
      testPostAccess("bam/json", "admin", testSample2) mustBe OK
      testPostAccess("bam/range", "", testSample1) mustBe OK
      testPostAccess("bam/range", "", testSample2) mustBe FORBIDDEN
      testPostAccess("bam/range", "admin", testSample1) mustBe OK
      testPostAccess("bam/range", "admin", testSample2) mustBe OK
      testPostAccess("bam/samtools", "", testSample1) mustBe OK
      testPostAccess("bam/samtools", "", testSample2) mustBe FORBIDDEN
      testPostAccess("bam/samtools", "admin", testSample1) mustBe OK
      testPostAccess("bam/samtools", "admin", testSample2) mustBe OK
    }

  }

  "`parseBamRequestFromPost`" should {

    "fail if there is no JSON body" in {
      val response = route(app, FakeAdminRequest(POST, "/bai")).get
      status(response) mustBe UNSUPPORTED_MEDIA_TYPE  // 415
      // Caught at the 'Action(parse.json)' level
    }

    "fail if there is no key 'sample' in the body" in {
      val body: JsValue = Json.parse(s"""{"yyyy": "xxx"}""")
      val response = route(app, FakeAdminRequest(POST, "/bai").withJsonBody(body)).get
      status(response) mustBe INTERNAL_SERVER_ERROR
      contentAsString(response) must startWith("No key found in request body")
    }

  }

  "`sampleNameToBamRequest`" should {

    "fail if the key is not known to the database (POST)" in {
      val body: JsValue = Json.parse(s"""{"sample": "xxx"}""")
      val response = route(app, FakeAuthorizedRequest(POST, "/bai").withJsonBody(body)).get
      status(response) mustBe INTERNAL_SERVER_ERROR
      contentAsString(response) must startWith("No corresponding BAM file")
    }

    "fail if the key exists but the bam file is not found on disk (POST)" in {
      val body: JsValue = Json.parse(s"""{"sample": "$notherekey"}""")
      val response = route(app, FakeAuthorizedRequest(POST, "/bai").withJsonBody(body)).get
      status(response) mustBe INTERNAL_SERVER_ERROR
      contentAsString(response) must startWith("This BAM file cannot be found")
    }

    "fail if the key is not known to the database (GET)" in {
      val response = route(app, FakeAuthorizedRequest(GET, s"/bai/xxx?token=$auth0Token")).get
      status(response) mustBe INTERNAL_SERVER_ERROR
      contentAsString(response) must startWith("No corresponding BAM file")
    }

    "fail if the key exists but the bam file is not found (GET)" in {
      val response = route(app, FakeAuthorizedRequest(GET, s"/bai/$notherekey?token=$auth0Token")).get
      status(response) mustBe INTERNAL_SERVER_ERROR
      contentAsString(response) must startWith("This BAM file cannot be found")
    }

  }
}