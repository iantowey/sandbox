package com.lucidlogic.horsetracker.model.dto;

import com.lucidlogic.horsetracker.config.Constants;
import com.lucidlogic.horsetracker.model.dto.MeetingDTO;
import com.parse.ParseClassName;
import com.parse.ParseObject;

import java.util.Date;
import java.util.List;

/**
 * Created by ian on 06/07/16.
 */
@ParseClassName(Constants.RACECARD)
public class RacecardDTO extends ParseObject{

    public Date getDate(){ return getDate(Constants.DATE); }
    public void  setDate(Date date){ put(Constants.DATE, date); }

    public List<MeetingDTO> getMeetings(){
        return getList(Constants.MEETINGS);
    }
    public void setMeetings(List<MeetingDTO> meetings){
        put(Constants.MEETINGS, meetings);
    }

}
