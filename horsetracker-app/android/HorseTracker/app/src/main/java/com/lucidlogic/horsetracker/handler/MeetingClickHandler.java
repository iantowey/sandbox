package com.lucidlogic.horsetracker.handler;

import android.view.View;

import com.lucidlogic.horsetracker.model.Meeting;
import com.lucidlogic.horsetracker.model.Race;

import timber.log.Timber;

/**
 * Created by ian on 12/07/16.
 */
public class MeetingClickHandler {

    public void onMeetingClick(View v, Meeting m){
        Timber.i("%s clicked" ,m.getTrack());
    }
}
