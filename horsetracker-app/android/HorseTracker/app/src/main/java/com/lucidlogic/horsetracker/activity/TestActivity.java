package com.lucidlogic.horsetracker.activity;

import android.support.v7.app.AppCompatActivity;
import android.os.Bundle;

import com.lucidlogic.horsetracker.R;
import com.lucidlogic.horsetracker.view.impl.RaceFragment;

public class TestActivity extends AppCompatActivity {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.content_main);
        this.getSupportFragmentManager().beginTransaction().replace(
                R.id.fragment_content,
                RaceFragment.newInstance("yhLeEP2v4m"),"RaceFragment").commit();
    }

}
