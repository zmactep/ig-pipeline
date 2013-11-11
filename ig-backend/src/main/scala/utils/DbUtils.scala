package utils

import java.sql._
import com.googlecode.protobuf.format.JsonFormat
import protocol.Command.ResponseBatch
import scala.collection.JavaConversions._
import java.io.IOException

/**
 * Created with IntelliJ IDEA.
 * User: Kos
 * Date: 29.10.13
 * Time: 9:19
 * To change this template use File | Settings | File Templates.
 */

object DbUtils {

  def saveWorkTask(work: Any, prep: PreparedStatement, conn: Connection): Int = {
    prep.setString(1, work.toString)
    prep.setString(2, "in progress")
    prep.setString(3, "in progress")
    prep.executeUpdate

    val statement = conn.createStatement()
    val resultSet = statement.executeQuery("SELECT LAST_INSERT_ID()")

    if (resultSet.next()) {
      resultSet.getInt(1)
    } else -1
  }

  def updateTask(id: Int, result: Any, prep: PreparedStatement, conn: Connection) = {
    var status: Boolean = true
    try {
      val builder = ResponseBatch.newBuilder()
      JsonFormat.merge(result.toString, builder)
      val responseBatch = builder.build()
      responseBatch.getResultList.foreach{result =>
        status = status && result.getStatus == "ok"
      }
    } catch {
      case e : IOException => status = false
    }

    prep.setString(1, result.toString)
    prep.setString(2, if (status) "ok" else "failed")
    prep.setString(3, id.toString)
    prep.executeUpdate
  }

  def findTask(id: Int, prep: PreparedStatement, conn: Connection) : String = {
    prep.setString(1, id.toString)
    val rs = prep.executeQuery()
    if (rs.next()) {
      rs.getString("result")
    } else {
      "Your job id not found"
    }
  }
}
