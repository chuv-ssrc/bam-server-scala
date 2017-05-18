package controllers.management

import javax.inject._

import auth.{AdminAction, AuthenticatedRequest}
import forms.SamplesForm
import models.Sample
import play.api.db.Database
import play.api.libs.json._
import play.api.mvc._

import scala.util.Try


@Singleton
class SamplesController @Inject()(db: Database) extends Controller {

  val AdminAction = new AdminAction(db)
  import SamplesController._

  /**
    * Add one or more samples to the database.
    * Expects a JSON body of the type
    * ```{ "samples": [{"name": "A", "filename": "/"}, {"name": "B", "filename": "/"}] }```
    */
  def addSamples() = AdminAction(parse.json) { implicit request =>

    val appId = request.user.appId
    val samples = samplesFromRequest

    val sampleCounts: Seq[Int] = samples.map(s => findSampleBySampleName(db, s.name, appId))
    if (sampleCounts.exists(_ > 0)) {
      val dupSample = samples(sampleCounts.indexWhere(_ > 0))
      InternalServerError(s"Cannot insert sample '${dupSample.name}' because it already exists. Nothing was inserted.")
    } else {
      db.withConnection { conn =>
        val unknowns: String = samples.map(_ => "(?,?)").mkString(",")
        val statement = conn.prepareStatement("INSERT INTO `samples`(`name`,`filename`) VALUES " + unknowns + " ;")
        samples.zipWithIndex.foreach { case (sample, i: Int) =>
          statement.setString(2 * i + 1, sample.name)
          statement.setString(2 * i + 2, sample.filename)
        }
        statement.execute()
        Ok(s"Inserted ${samples.size} sample(s)")
      }
    }
  }

  /**
    * Delete one or more samples from the database.
    * Expects a JSON body of the type
    * ```{"samples": ["A","B"]}```
    */
  def deleteSamples() = AdminAction(parse.json) { implicit request =>

    val appId = request.user.appId
    val sampleNames = sampleNamesFromRequest

    val sampleCounts: Seq[Int] = sampleNames.map(sname => findSampleBySampleName(db, sname, appId))
    if (sampleCounts.contains(0)) {
      val unknownSampleName = sampleNames(sampleCounts.find(_ == 0).get)
      InternalServerError(s"Cannot delete sample '$unknownSampleName' because it does not exist. Nothing was deleted.")
    } else {
      db.withConnection { conn =>
        val unknowns: String = sampleNames.map(_ => "?").mkString(",")
        val statement = conn.prepareStatement("DELETE FROM `samples` WHERE `name` IN (" + unknowns + ") ;")
        sampleNames.zipWithIndex.foreach { case (name, i: Int) => statement.setString(i + 1, name) }
        statement.execute()
        Ok(s"Deleted ${sampleNames.size} sample(s)")
      }
    }
  }

}


object SamplesController {

  def samplesFromRequest(implicit request: AuthenticatedRequest[JsValue]): Seq[Sample] = {
    val samplesJs: JsArray = (request.body \ "samples").asOpt[JsArray] getOrElse {
      throw new IllegalArgumentException("Could not cast samples array from request body to JsArray")
    }
    val samples: Seq[Sample] = samplesJs.value map { sampleJs =>
      Try (SamplesForm.fromJson(sampleJs)) getOrElse {
        throw new IllegalArgumentException("Could not cast sample from request body to Samples model")
      }
    }
    samples
  }

  def sampleNamesFromRequest(implicit request: AuthenticatedRequest[JsValue]): Seq[String] = {
    val sampleNamesJs: JsArray = (request.body \ "samples").asOpt[JsArray] getOrElse {
      throw new IllegalArgumentException("Could not cast sample names array from request body to JsArray")
    }
    val sampleNames: Seq[String] = sampleNamesJs.value map { sampleNameJs =>
      Try (sampleNameJs.as[String]) getOrElse {
        throw new IllegalArgumentException("Could not cast sample names from request body to String")
      }
    }
    sampleNames
  }

  /**
    * Return the number of samples found with this *sampleName* and *appId* in the database.
    */
  def findSampleBySampleName(db: Database, sampleName: String, appId: Int): Int = {
    db.withConnection { conn =>
      val statement = conn.prepareStatement("SELECT `id` FROM `samples` WHERE name = ? ;")
      statement.setString(1, sampleName)
      //statement.setInt(2, appId)
      val result = statement.executeQuery()
      var count = 0
      while (result.next()) {
        count += 1
      }
      count
    }
  }


}
