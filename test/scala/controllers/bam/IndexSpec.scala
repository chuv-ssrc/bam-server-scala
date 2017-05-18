package scala.controllers.bam

import org.scalatestplus.play._
import play.api.libs.json._
import play.api.test.Helpers._
import play.api.test._

import scala.setup.WithToken


/**
  * Test IndexController.
  */
class IndexSpec extends PlaySpec with OneAppPerSuite with WithToken {

  "IndexController" should {

    "provide the BAM index if everything is right (POST)" in {
      val body: JsValue = Json.parse(s"""{"sample": "$testSample1"}""")
      val response = route(app, FakeAuthorizedRequest(POST, "/bai").withJsonBody(body)).get
      status(response) mustBe OK
      contentType(response).get mustBe BINARY  // application/octet-stream
    }

    "provide the BAM index if everything is right (GET)" in {
      val response = route(app, FakeAuthorizedRequest(GET, s"/bai/$testSample1")).get
      status(response) mustBe OK
      contentType(response).get mustBe BINARY
    }

    "provide the BAM index if everything is right, with token in URL (GET)" in {
      val response = route(app, FakeRequest(GET, s"/bai/$testSample1?token=$auth0Token")).get
      status(response) mustBe OK
      contentType(response).get mustBe BINARY
    }

  }

}
