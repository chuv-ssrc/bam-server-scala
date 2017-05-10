package scala.controllers.management

import controllers.management.AppsController
import org.scalatestplus.play._
import play.api.libs.json._
import play.api.test.Helpers._
import utils.DbUtils._

import scala.setup.{WithDatabase, WithToken}


/**
  * Test AppsController.
  */
class AppsSpec extends PlaySpec with OneAppPerSuite with WithToken with WithDatabase {

  val db = dbContext.db

  "`findAppByIss`" should {

    "return 0 if the app does not exist" in {
      AppsController.findAppByIss(db, "unknown") must equal(0)
    }

    "return 1 if the app already exists" in {
      AppsController.findAppByIss(db, "test") must equal(1)
    }

  }

  "Add an app" should {

    "add one rows to the app table" in {
      val body = Json.parse("""
        { "iss": "testApp", "keyFile": "/", "description": "none" }
      """)
      val db = dbContext.db
      val count0 = countRows(db, "apps")
      val response = route(app, FakeAdminRequest(PUT, "/apps").withJsonBody(body)).get
      status(response) mustBe OK
      val count = countRows(db, "apps")
      count must equal(count0 + 1)
    }

    "fail if the app already exists" in {
      val body = Json.parse("""
        { "iss": "testApp", "keyFile": "/", "description": "none" }
      """)
      val db = dbContext.db
      val count0 = countRows(db, "apps")
      val response = route(app, FakeAdminRequest(PUT, "/apps").withJsonBody(body)).get
      status(response) mustBe INTERNAL_SERVER_ERROR
      val count = countRows(db, "apps")
      count must equal(count0)
    }

  }

  "Delete an app" should {

    "remove one rows from the apps table" in {
      val body = Json.parse("""
        { "iss": "testApp" }
      """)
      val db = dbContext.db
      val count0 = countRows(db, "apps")
      val response = route(app, FakeAdminRequest(DELETE, "/apps").withJsonBody(body)).get
      status(response) mustBe OK
      val count = countRows(db, "apps")
      count must equal(count0 - 1)
    }

    "fail if the app does not exist" in {
      val body = Json.parse("""
        { "iss": "testApp" }
      """)
      val db = dbContext.db
      val count0 = countRows(db, "apps")
      val response = route(app, FakeAdminRequest(DELETE, "/apps").withJsonBody(body)).get
      status(response) mustBe INTERNAL_SERVER_ERROR
      val count = countRows(db, "apps")
      count must equal(count0)
    }

  }

}
