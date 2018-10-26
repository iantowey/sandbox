package com.lucidlogic.horsetracker.adapter;

import android.support.v4.app.Fragment;
import android.support.v4.view.PagerAdapter;
import android.view.View;

import java.util.List;

/**
 * Created by ian on 28/07/16.
 */
public class RacecardViewPagerAdapter extends PagerAdapter {

    private List<Fragment> fragmentList;

    public RacecardViewPagerAdapter(List<Fragment> fragmentList) {
        this.fragmentList = fragmentList;
    }


    @Override
    public int getCount() {
        return fragmentList.size();
    }

    @Override
    public boolean isViewFromObject(View view, Object object) {
        return false;
    }

}
