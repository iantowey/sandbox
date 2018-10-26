package com.lucidlogic.horsetracker;

import android.app.Application;

import com.lucidlogic.horsetracker.config.AppConfig;
import com.lucidlogic.horsetracker.service.ParseService;
import com.parse.ParseException;
import com.parse.ParseUser;

/**
 * Created by ian on 06/07/16.
 */
public class App extends Application {

    @Override
    public void onCreate() {
        super.onCreate();
        AppConfig.init(this);

        if(ParseUser.getCurrentUser() == null){
            try {
                ParseUser.logIn("ian", "ian");
            } catch (ParseException e) {
                e.printStackTrace();
            }
        }
        ParseService.initUser();

    }
}
